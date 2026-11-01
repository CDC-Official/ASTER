// RKF78 专项验证：二体解析解、误差拒绝、失败保留状态及正反向终点。
// 默认仅需临时合成 EOP；传入 --jpl <文件> 可额外验证 EBCRS/TDB 接口。
#include "core_propagation.hpp"
#include <array>
#include <chrono>
#include <filesystem>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace {
namespace fs = std::filesystem;
using State = std::array<double, 6>;
using bastion::Overwatch;
int checks = 0;

void require(bool condition, const char* message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

// 独立二体参考：用长双精度解开普勒方程，再直接旋转轨道平面状态；不调用项目转换/传播函数。
State analytic(double elapsed)
{
    const long double a = 2.2L, e = .15L, inc = .4L, node = .3L, peri = .2L;
    const long double M = .5L + static_cast<long double>(elapsed) / std::sqrt(a*a*a);
    long double E = M;
    for (int k = 0; k < 30; ++k) {
        const long double correction = (E-e*std::sin(E)-M)/(1-e*std::cos(E));
        E -= correction;
        if (std::abs(correction) < 1e-18L) break;
    }
    const long double cN=std::cos(node), sN=std::sin(node), cw=std::cos(peri), sw=std::sin(peri);
    const long double ci=std::cos(inc), si=std::sin(inc), q=std::sqrt(1-e*e);
    const long double P[3]={cN*cw-sN*sw*ci, sN*cw+cN*sw*ci, sw*si};
    const long double Q[3]={-cN*sw-sN*cw*ci, -sN*sw+cN*cw*ci, cw*si};
    const long double rX=a*(std::cos(E)-e), rY=a*q*std::sin(E);
    const long double vX=-std::sin(E)/(std::sqrt(a)*(1-e*std::cos(E)));
    const long double vY=q*std::cos(E)/(std::sqrt(a)*(1-e*std::cos(E)));
    State result{};
    for (int k=0;k<3;++k) {
        result[k]=static_cast<double>(P[k]*rX+Q[k]*rY);
        result[k+3]=static_cast<double>(P[k]*vX+Q[k]*vY);
    }
    return result;
}

double distance(const State& x, const State& y, int offset)
{
    double sum=0;
    for(int k=offset;k<offset+3;++k) sum+=(x[k]-y[k])*(x[k]-y[k]);
    return std::sqrt(sum);
}

void check_trace(const Overwatch& tracer, const std::string& function)
{
    std::ostringstream trace,events;
    tracer.ReportTrace(trace); tracer.ReportEvents(events);
    require(trace.str().find(function)!=std::string::npos,"Failure lost integrator call stack");
    require(events.str().find("[ERROR]")!=std::string::npos,"Failure missing event log error");
}

// GCRS 与 E-BCRS 的时标不同，但共用同一组 RKF 参数与误差／失败检查。
template<class Job>
auto& rkfSetting(Job& job)
{
    return job.opn_setting.rkf;
}

// 所有单步接口共用这些输入/失败检查。失败后原数组与 h_return 均检查，不仅检查错误码。
template<class Job, class Step>
void check_steps(Job job, Step step, const State& initial, const char* name)
{
    auto fails=[&](Job invalid, State x, double epoch, double h, const char* message) {
        State original=x; double returned=123;
        Overwatch tracer(__FILE__,"expected_failure",__LINE__);
        step(tracer,invalid,epoch,x.data(),h,returned);
        require(tracer.Code()!=ERR_NOERROR,message);
        require(returned==0,"Failed step returned nonzero accepted time");
        for(int k=0;k<6;++k)
            require(x[k]==original[k] || (std::isnan(x[k]) && std::isnan(original[k])),"Failed step changed initial state");
        check_trace(tracer,name);
    };
    // 强制触发试算上限。对 NSG2 增加 J2 扰动，否则二体下仅平近点角匀速，误差为零。
    Job strict=job; strict.opn_model.use_earth_NSG=true;
    strict.opn_model.earth_zonal_order=2; strict.opn_model.earth_tesseral_degree=0;
    rkfSetting(strict).maxTE_norm=1e-30; rkfSetting(strict).iterCount_max=1;
    fails(strict,initial,60000.25,.5,"Exhausted step accepted as success");
    fails(strict,initial,60000.25,-.5,"Negative step bypassed error control");

    // 拒绝较大步长后，接受的结果必须等于从头用最终步长重算的结果。
    Overwatch retry_trace(__FILE__,"retry_success",__LINE__);
    rkfSetting(strict).maxTE_norm=1e-12; rkfSetting(strict).iterCount_max=30;
    State retried=initial; double accepted=0;
    step(retry_trace,strict,60000.25,retried.data(),2.0,accepted);
    require(retry_trace.Code()==ERR_NOERROR,"Halving failed to recover a usable step");
    require(accepted>0 && accepted<2.0,"Test did not exercise a rejected step");
    State direct=initial; double direct_h=0;
    Overwatch direct_trace(__FILE__,"direct_success",__LINE__);
    step(direct_trace,strict,60000.25,direct.data(),accepted,direct_h);
    require(direct_trace.Code()==ERR_NOERROR && direct_h==accepted,"Accepted step not reproducible");
    require(retried==direct,"Accepted step used derivatives from a different step size");

    for(double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        Job invalid=job;rkfSetting(invalid).maxTE_norm=bad;
        fails(invalid,initial,60000.25,.1,"Invalid tolerance accepted");
    }
    for(int bad:{0,-1}) {
        Job invalid=job;rkfSetting(invalid).iterCount_max=bad;
        fails(invalid,initial,60000.25,.1,"Invalid trial limit accepted");
    }
    for(double h:{0.,1e-20,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        fails(job,initial,60000.25,h,"Invalid/nonadvancing step accepted");
    fails(job,initial,std::numeric_limits<double>::quiet_NaN(),.1,"NaN epoch accepted");
    State bad=initial;bad[0]=std::numeric_limits<double>::quiet_NaN();
    fails(job,bad,60000.25,.1,"NaN state accepted");
    Job missing=job;missing.opn_EOP=EOP{};
    fails(missing,initial,60000.25,.1,"Derivative failure not propagated");
    if(std::string(name)!="RKF78_NSG2") {
        // 零半径使二体加速度非有限；应返回错误，而不是接受 NaN 状态。
        bad=initial;bad[0]=bad[1]=bad[2]=0;
        fails(job,bad,60000.25,.1,"Degenerate Cartesian state accepted");
    }
    double returned=123;Overwatch null_trace(__FILE__,"null_state",__LINE__);
    step(null_trace,job,60000.25,nullptr,.1,returned);
    require(null_trace.Code()!=ERR_NOERROR && returned==0,"Null state accepted");
}

// 同一接受步的误差必须可重现；失败返回 NaN，不能残留上次成功的误差。
template<class Job, class Step>
void check_monitored(Job job, Step step, const State& initial)
{
    job.opn_setting.display_progress=0;
    rkfSetting(job).maxTE_norm=1e-12;
    rkfSetting(job).iterCount_max=30;
    job.opn_model.use_earth_NSG=true;
    job.opn_model.earth_zonal_order=2;
    job.opn_model.earth_tesseral_degree=0;
    State state=initial;double h=0,error=-1;
    Overwatch trace(__FILE__,"accepted_error",__LINE__);
    step(trace,job,60000.25,state.data(),2.,h,error);
    require(trace.Code()==ERR_NOERROR && h>0 && h<2.,"Monitored trial did not recover after rejection");
    require(std::isfinite(error) && error>=0 && error<=rkfSetting(job).maxTE_norm,"Invalid accepted error");
    State direct=initial;double dh=0,de=-1;
    Overwatch direct_trace(__FILE__,"direct_error",__LINE__);
    step(direct_trace,job,60000.25,direct.data(),h,dh,de);
    require(direct_trace.Code()==ERR_NOERROR && direct==state && dh==h && de==error,"Error/state changed on direct accepted step");
    rkfSetting(job).maxTE_norm=1e-30;rkfSetting(job).iterCount_max=1;
    direct=initial;error=123;
    Overwatch failure(__FILE__,"failed_error",__LINE__);
    step(failure,job,60000.25,direct.data(),2.,dh,error);
    require(failure.Code()!=ERR_NOERROR && dh==0 && std::isnan(error) && direct==initial,"Failed step returned stale error/state");
}

// 定步长控制独立于容差：保留原 RKF 公式，只改变超限后的处理。
template<class Job, class Step>
void check_fixed_steps(Job job, Step step, const State& initial)
{
    job.opn_model.use_earth_NSG=true;
    job.opn_model.earth_zonal_order=2;
    job.opn_model.earth_tesseral_degree=0;
    rkfSetting(job).fixed_step=true;
    job.opn_setting.display_progress=0;
    rkfSetting(job).maxTE_norm=1e-30;
    rkfSetting(job).iterCount_max=1;
    for(double input_h:{.5,-.5}) {
        State state=initial;double h=0,error=-1;
        Overwatch trace(__FILE__,"fixed_step",__LINE__);
        step(trace,job,60000.25,state.data(),input_h,h,error);
        require(trace.Code()==ERR_NOERROR && h==input_h,"Fixed mode changed the requested step");
        require(std::isfinite(error) && error>rkfSetting(job).maxTE_norm,"Fixed-mode test did not exceed tolerance");
        std::ostringstream events;trace.ReportEvents(events);
        require(events.str().find("[WARNING]")!=std::string::npos,"Fixed mode lost over-tolerance warning");
        require(events.str().find("fixed step retained")!=std::string::npos,"Fixed warning does not explain retained step");
        // 与宽松容差下的一次原公式计算逐分量相同，排除改换公式或用旧导数推进。
        Job loose=job;rkfSetting(loose).fixed_step=false;rkfSetting(loose).maxTE_norm=1.;
        State reference=initial;double rh=0,re=-1;
        Overwatch comparison(__FILE__,"fixed_formula_reference",__LINE__);
        step(comparison,loose,60000.25,reference.data(),input_h,rh,re);
        require(comparison.Code()==ERR_NOERROR && rh==h && reference==state && re==error,"Fixed mode changed the RKF formula or error estimate");
    }
    Job missing=job;missing.opn_EOP=EOP{};
    State failed=initial;double h=1,error=0;
    Overwatch failure(__FILE__,"fixed_derivative_failure",__LINE__);
    step(failure,missing,60000.25,failed.data(),.5,h,error);
    require(failure.Code()!=ERR_NOERROR && failed==initial && h==0 && std::isnan(error),"Fixed mode accepted a force-model failure");
}

template<class Job, class Propagate, class Orbit>
void check_propagation(Job job, Propagate propagate, bool tdb, const char* name)
{
    // 归一化 μ=1 的椭圆轨道；TT/TDB 内的时间长度用于解析解，GCRS 输出历元转换为 UTC。
    double max_position=0,max_velocity=0,max_endpoint=0;
    for(double duration:{7.13,-7.13,.11,-.11}) {
        job.opn_setting.duration_norm=duration;
        job.opn_setting.h0_norm=.3;
        rkfSetting(job).maxTE_norm=1e-12;
        rkfSetting(job).iterCount_max=20;
        std::vector<double> epochs;std::vector<Orbit> states;
        Overwatch tracer(__FILE__,"analytic_propagation",__LINE__);
        propagate(tracer,job,epochs,states);
        require(tracer.Code()==ERR_NOERROR,"Analytic propagation failed");
        require(epochs.size()==states.size() && epochs.size()>=2,"Incomplete propagated output");
        State last{};
        std::array<double, 6> exported_values_1; states.back().ToArray(exported_values_1.data());
        if constexpr(std::is_same_v<Orbit,Cartesian>) last=exported_values_1;
        else { states.back().ToCartesian(tracer,last.data()); require(tracer.Code()==ERR_NOERROR,"NSG2 final conversion failed"); }
        const State reference=analytic(duration);
        max_position=std::max(max_position,distance(last,reference,0)*ae);
        max_velocity=std::max(max_velocity,distance(last,reference,3)*NORMALIZED_VELOCITY_TO_MPS);
        std::cout<<"CASE "<<name<<" duration="<<duration<<" position_m="<<distance(last,reference,0)*ae<<" velocity_mps="<<distance(last,reference,3)*NORMALIZED_VELOCITY_TO_MPS<<'\n';
        require(distance(last,reference,0)<2e-9,"Two-body position disagrees with analytic reference");
        require(distance(last,reference,3)<2e-9,"Two-body velocity disagrees with analytic reference");
        const double epoch_elapsed=(epochs.back()-epochs.front())*86400.;
        const double endpoint_error=std::abs(epoch_elapsed-duration*NORMALIZED_TIME_TO_SEC);
        max_endpoint=std::max(max_endpoint,endpoint_error);
        // UTC 包装函数先合并为完整 JD，会引入数十微秒舍入；不将该历史误差归因于积分器。
        const double expected_epoch=tdb ? 60000.25+duration*NORMALIZED_TIME_TO_DAY
            : TT2UTC(UTC2TT(60000.25)+duration*NORMALIZED_TIME_TO_DAY);
        require(std::abs(epochs.back()-expected_epoch)*86400.<1e-6,"Endpoint clipping failed");
        require(endpoint_error<(tdb?1e-6:1e-4),"Recorded epoch exceeds existing time-conversion precision");
        for(size_t k=1;k<epochs.size();++k)
            require((epochs[k]-epochs[k-1])*duration>0,"Recorded epochs moved in wrong direction");
        // 整段预报的步长设置只取大小，负步长也应得到相同的有效结果。
        job.opn_setting.h0_norm=-.3;
        std::vector<double> negative_epochs;std::vector<Orbit> negative_states;
        Overwatch negative_trace(__FILE__,"negative_step_magnitude",__LINE__);
        propagate(negative_trace,job,negative_epochs,negative_states);
        require(negative_trace.Code()==ERR_NOERROR && epochs==negative_epochs,"Step sign changed propagation direction");
        std::array<double, 6> exported_values_2; negative_states.back().ToArray(exported_values_2.data());
        std::array<double, 6> exported_values_3; states.back().ToArray(exported_values_3.data());
        require(exported_values_3==exported_values_2,"Step magnitude policy changed results");
    }
    // 零时长仅记录初态；坏设置不应添加初态，单步失败只保留已接受的初态。
    job.opn_setting.duration_norm=0;job.opn_setting.h0_norm=0;
    std::vector<double> epochs;std::vector<Orbit> states;
    Overwatch zero_trace(__FILE__,"zero_duration",__LINE__);
    propagate(zero_trace,job,epochs,states);
    require(zero_trace.Code()==ERR_NOERROR && epochs.size()==1 && states.size()==1,"Zero duration did not retain only initial state");
    job.opn_setting.duration_norm=1;
    Overwatch invalid_trace(__FILE__,"zero_step",__LINE__);
    propagate(invalid_trace,job,epochs,states);
    require(invalid_trace.Code()!=ERR_NOERROR && epochs.size()==1,"Invalid settings modified existing output");
    job.opn_setting.h0_norm=.5;rkfSetting(job).maxTE_norm=1e-30;rkfSetting(job).iterCount_max=1;
    job.opn_model.use_earth_NSG=true;job.opn_model.earth_zonal_order=2;job.opn_model.earth_tesseral_degree=0;
    Overwatch failed_trace(__FILE__,"failed_outer_step",__LINE__);
    propagate(failed_trace,job,epochs,states);
    require(failed_trace.Code()!=ERR_NOERROR && epochs.size()==2 && states.size()==2,"Rejected state appended to ephemeris");
    // 跨出 EOP 覆盖区间，在已经完成若干有效步之后触发右函数错误。
    // 输出应停在最后有效记录，不能清空旧点，也不能追加失败步骤。
    Job partial=job;
    partial.opn_setting.duration_norm=1;
    partial.opn_setting.h0_norm=.1;
    rkfSetting(partial).maxTE_norm=1e-12;
    rkfSetting(partial).iterCount_max=20;
    if constexpr(std::is_same_v<Job,OrbitPropagationJob>) partial.opn_setting.t0_MJD_UTC=60001.998;
    else partial.opn_setting.t0_MJD_TDB=60001.998;
    std::vector<double> partial_epochs;std::vector<Orbit> partial_states;
    Overwatch partial_trace(__FILE__,"partial_ephemeris",__LINE__);
    propagate(partial_trace,partial,partial_epochs,partial_states);
    require(partial_trace.Code()!=ERR_NOERROR,"Propagation outside EOP coverage succeeded");
    require(partial_states.size()>1 && partial_states.size()==partial_epochs.size(),"Lost accepted records on a later failure");
    require(partial_epochs.back()<60002.,"Failed step appended outside EOP coverage");
    for(const auto& state:partial_states) {
        std::array<double, 6> exported_values_4; state.ToArray(exported_values_4.data());
        for(double value:exported_values_4) require(std::isfinite(value),"Partial ephemeris contains invalid state");
    }
    std::cout << "METRIC " << name << " position_m=" << max_position << " velocity_mps=" << max_velocity
              << " endpoint_s=" << max_endpoint << " time_system=" << (tdb?"TDB":"TT/UTC") << '\n';
}
}

int main(int argc,char** argv) try {
    const fs::path root=fs::temp_directory_path()/("cdc-rkf78-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(root);
    struct Cleanup {fs::path path;~Cleanup(){std::error_code ec;fs::remove_all(path,ec);}} cleanup{root};
    // 覆盖正反向积分时段的连续 EOP，仅供本测试，不是生产地球定向数据。
    const auto eop_file=root/"eop.txt";
    std::ofstream eop(eop_file);
    eop << "BEGIN OBSERVED\n"
        << "2023 2 24 59999 0.1 0.2 0.3 0 0 0 0 0 37\n"
        << "2023 2 25 60000 0.1 0.2 0.3 0 0 0 0 0 37\n"
        << "2023 2 26 60001 0.1 0.2 0.3 0 0 0 0 0 37\n"
        << "2023 2 27 60002 0.1 0.2 0.3 0 0 0 0 0 37\nEND OBSERVED\n";
    eop.close();require(bool(eop),"Cannot write EOP fixture");
    OrbitPropagationJob job;
    Overwatch setup(__FILE__,"setup",__LINE__);
    std::string eop_path=eop_file.string();
    job.opn_EOP.Read_EOP(setup,eop_path);
    job.opn_setting.t0_MJD_UTC=60000.25;
    const State initial=analytic(0);
    job.initCart_norm.Set(initial[0],initial[1],initial[2],initial[3],initial[4],initial[5],true);
    job.initCart_norm.ToKeplerian(setup,job.initKep_norm);
    NonSingular2 nsg;job.initCart_norm.ToNonSingular2(setup,nsg);
    require(setup.Code()==ERR_NOERROR,"Cannot initialize two-body test orbit");
    std::cout<<std::setprecision(17);
    if(argc==2 && std::string(argv[1])=="--probe") {
        rkfSetting(job).maxTE_norm=1e-30;rkfSetting(job).iterCount_max=1;
        State x=initial;double h=123;Overwatch probe(__FILE__,"probe",__LINE__);
        RKF78_Cart(probe,job,60000.25,x.data(),.5,h);
        std::cout<<"PROBE error="<<probe.Code()<<" returned_h="<<h<<" unchanged="<<(x==initial)<<'\n';
        return 0;
    }
    // 单独验证 NSG2 跨圈转换，覆盖广义偏近点经度与平经度位于周界两侧的情况。
    const bool wrap_probe=argc==2 && std::string(argv[1])=="--wrap-probe";
    for(double elapsed:{0.,-7.13}) {
        std::array<double, 6> exported_values_5; nsg.ToArray(exported_values_5.data());
        State base_elements=exported_values_5;base_elements[5]+=elapsed/std::sqrt(2.2*2.2*2.2);
        const State expected=analytic(elapsed);
        for(int turns:{-3,-1,0,1,3}) {
            State elements=base_elements;elements[5]+=turns*2*pi;
            Overwatch wrap_trace(__FILE__,"nsg2_wrap",__LINE__);
            NonSingular2 shifted(wrap_trace,elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],true);
            State restored{};shifted.ToCartesian(wrap_trace,restored.data());
            const double position_error=distance(restored,expected,0)*ae;
            if(wrap_probe) std::cout<<"WRAP elapsed="<<elapsed<<" turns="<<turns<<" position_m="<<position_error<<'\n';
            else {
                require(wrap_trace.Code()==ERR_NOERROR,"NSG2 wrap conversion failed");
                require(position_error<1e-5,"NSG2 angle wrap changed the physical state");
            }
        }
    }
    if(wrap_probe) return 0;
    // 显式选取原六参数接口，继续覆盖既有调用方式；下面另验带误差输出的重载。
    using CartesianStep = void (*)(Overwatch&,const OrbitPropagationJob&,double,double*,double,double&);
    check_steps(job,static_cast<CartesianStep>(RKF78_Cart),initial,"RKF78_Cart");
    {
        OrbitPropagationJob monitored=job;
        monitored.opn_model.use_earth_NSG=true;
        monitored.opn_model.earth_zonal_order=2;
        monitored.opn_model.earth_tesseral_degree=0;
        monitored.opn_setting.display_progress=0;
        rkfSetting(monitored).maxTE_norm=1e-12;
        rkfSetting(monitored).iterCount_max=30;
        State retried=initial; double accepted=0, estimate=-1;
        Overwatch trace(__FILE__,"monitored_step",__LINE__);
        RKF78_Cart(trace,monitored,60000.25,retried.data(),2.,accepted,estimate);
        require(trace.Code()==ERR_NOERROR && accepted>0 && accepted<2.,"Monitored step did not recover after rejection");
        require(std::isfinite(estimate) && estimate>=0 && estimate<=rkfSetting(monitored).maxTE_norm,"Invalid accepted-step error estimate");
        State direct=initial; double direct_h=0,direct_estimate=-1;
        Overwatch direct_trace(__FILE__,"monitored_direct",__LINE__);
        RKF78_Cart(direct_trace,monitored,60000.25,direct.data(),accepted,direct_h,direct_estimate);
        require(direct_trace.Code()==ERR_NOERROR && direct==retried && direct_estimate==estimate,"Displayed estimate does not belong to accepted step");
        rkfSetting(monitored).maxTE_norm=1e-30;rkfSetting(monitored).iterCount_max=1;
        State failed=initial;double failed_h=1,failed_estimate=0;
        Overwatch failed_trace(__FILE__,"monitored_failure",__LINE__);
        RKF78_Cart(failed_trace,monitored,60000.25,failed.data(),.5,failed_h,failed_estimate);
        require(failed_trace.Code()!=ERR_NOERROR && failed==initial && failed_h==0 && std::isnan(failed_estimate),"Failed step leaked a successful error estimate");
    }
    std::array<double, 6> exported_values_6; nsg.ToArray(exported_values_6.data());
    check_steps(job,static_cast<CartesianStep>(RKF78_NSG2),exported_values_6,"RKF78_NSG2");
    using MonitoredEarth = void (*)(Overwatch&,const OrbitPropagationJob&,double,double*,double,double&,double&);
    std::array<double, 6> exported_values_7; nsg.ToArray(exported_values_7.data());
    check_monitored(job,static_cast<MonitoredEarth>(RKF78_NSG2),exported_values_7);
    std::array<double, 6> exported_values_8; nsg.ToArray(exported_values_8.data());
    check_fixed_steps(job,static_cast<MonitoredEarth>(RKF78_NSG2),exported_values_8);
    check_fixed_steps(job,static_cast<MonitoredEarth>(RKF78_Cart),initial);
    check_propagation<OrbitPropagationJob,decltype(&OrbitPropagation_RKF78_Cart),Cartesian>(job,OrbitPropagation_RKF78_Cart,false,"GCRS-Cart");
    check_propagation<OrbitPropagationJob,decltype(&OrbitPropagation_RKF78_NSG2),NonSingular2>(job,OrbitPropagation_RKF78_NSG2,false,"GCRS-NSG2");
    if(argc==3 && std::string(argv[1])=="--jpl") {
        require(Read_JPLEPH(argv[2])==1,"Cannot load optional JPL ephemeris");
        OrbitPropagationJob_EBCRS ebcrs;
        ebcrs.opn_EOP=job.opn_EOP;ebcrs.initCart_norm=job.initCart_norm;ebcrs.initKep_norm=job.initKep_norm;
        ebcrs.opn_setting.t0_MJD_TDB=60000.25;
        using EbcrsStep = void (*)(Overwatch&,const OrbitPropagationJob_EBCRS&,double,double*,double,double&);
        check_steps(ebcrs,static_cast<EbcrsStep>(RKF78_Cart_EBCRS),initial,"RKF78_Cart_EBCRS");
        using MonitoredEbcrs = void (*)(Overwatch&,const OrbitPropagationJob_EBCRS&,double,double*,double,double&,double&);
        check_monitored(ebcrs,static_cast<MonitoredEbcrs>(RKF78_Cart_EBCRS),initial);
        check_fixed_steps(ebcrs,static_cast<MonitoredEbcrs>(RKF78_Cart_EBCRS),initial);
        check_propagation<OrbitPropagationJob_EBCRS,decltype(&OrbitPropagation_RKF78_Cart_EBCRS),Cartesian>(ebcrs,OrbitPropagation_RKF78_Cart_EBCRS,true,"EBCRS-Cart");
    }
    std::cout<<"RKF78 checks passed: "<<checks<<'\n';
    return 0;
} catch(const std::exception& error) {
    std::cerr<<"RKF78 test failed: "<<error.what()<<'\n';return 1;
}
