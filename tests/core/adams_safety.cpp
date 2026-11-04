#include <array>
// Adams 专项：导师原核逐步对照、解析解、系数、历史重置与失败原子性。
#include "core_propagation.hpp"
#include "accel_complete.hpp"
#include "nao_adams_reference.hpp"
#include <chrono>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
using State = std::array<double,6>;
namespace {
int checks=0;
void require(bool condition,const char* message) {
    ++checks; if(!condition) throw std::runtime_error(message);
}
bool equalHistory(const AdamsHistory& a,const AdamsHistory& b) {
    return a.steps==b.steps && a.filled==b.filled && a.use_correction==b.use_correction &&
        a.step_norm==b.step_norm && a.epoch_MJD==b.epoch_MJD && a.elapsed_norm==b.elapsed_norm &&
        a.state_norm==b.state_norm && a.derivatives==b.derivatives &&
        a.predictor==b.predictor && a.corrector==b.corrector;
}
// 长双精度的独立二体解析状态，不调用本库传播或状态转换函数。
State analytic(double elapsed) {
    const long double a=2.2L,e=.15L,inc=.4L,node=.3L,peri=.2L;
    long double M=.5L+elapsed/std::sqrt(a*a*a),E=M;
    for(int k=0;k<40;++k) {long double d=(E-e*std::sin(E)-M)/(1-e*std::cos(E));E-=d;if(std::abs(d)<1e-18L)break;}
    const long double P[3]={std::cos(node)*std::cos(peri)-std::sin(node)*std::sin(peri)*std::cos(inc),
        std::sin(node)*std::cos(peri)+std::cos(node)*std::sin(peri)*std::cos(inc),std::sin(peri)*std::sin(inc)};
    const long double Q[3]={-std::cos(node)*std::sin(peri)-std::sin(node)*std::cos(peri)*std::cos(inc),
        -std::sin(node)*std::sin(peri)+std::cos(node)*std::cos(peri)*std::cos(inc),std::cos(peri)*std::sin(inc)};
    long double r[2]={a*(std::cos(E)-e),a*std::sqrt(1-e*e)*std::sin(E)};
    long double v[2]={-std::sin(E)/(std::sqrt(a)*(1-e*std::cos(E))),std::sqrt(1-e*e)*std::cos(E)/(std::sqrt(a)*(1-e*std::cos(E)))};
    State x{};for(int k=0;k<3;++k){x[k]=P[k]*r[0]+Q[k]*r[1];x[k+3]=P[k]*v[0]+Q[k]*v[1];}return x;
}
double separation(const State& a,const State& b,int offset) {
    double s=0;for(int i=offset;i<offset+3;++i)s+=(a[i]-b[i])*(a[i]-b[i]);return std::sqrt(s);
}

void check_reference(OrbitPropagationJob& job,const State& initial) {
    job.opn_model.use_earth_NSG=true;
    job.opn_model.earth_zonal_order=2; job.opn_model.earth_tesseral_degree=0;
    const double t0=UTC2TT(job.opn_setting.t0_MJD_UTC);
    for(int k=2;k<=9;++k) for(bool correction:{false,true}) for(double h:{.01,-.01}) {
        job.opn_setting.adams={k,correction};
        State x=initial;AdamsHistory history;
        double elapsed=0,accepted,error;
        bastion::Overwatch trace(__FILE__,"reference",__LINE__);
        // 用生产入口起步，并确认缓存的是接受状态的真实末端导数。
        for(int i=0;i<k;++i) {
            Adams_NSG2(trace,job,t0,elapsed,x.data(),h,history,accepted,error);
            require(trace.Code()==ERR_NOERROR && accepted==h && std::isfinite(error),"Startup failed");
            elapsed+=h;
            State f{};Acceleration_NSG2(trace,job,t0+elapsed*constants::Conversion::NORMALIZED_TIME_TO_DAY,6,x.data(),f.data());
            require(history.derivatives[i]==f,"Startup derivative not evaluated at accepted state");
        }
        // 用同一批历史节点，直接执行导师原函数与新函数；不借两个力模型的差异解释结果。
        nao_reference::OrbitPropagationNModel model{job};
        nao_reference::d_matrix y(k,std::vector<double>(6));
        for(int j=0;j<k;++j)std::copy(history.derivatives[j].begin(),history.derivatives[j].end(),y[j].begin());
        State reference=x;double reference_time=elapsed;
        for(int i=0;i<20;++i) {
            nao_reference::Step(trace,6,k,h,t0,model,reference_time,reference.data(),y,correction);
            Adams_NSG2(trace,job,t0,elapsed,x.data(),h,history,accepted,error);
            require(trace.Code()==ERR_NOERROR && accepted==h && std::isnan(error),"Normal Adams step failed");
            require(x==reference,"Port differs from original NAO kernel");
            for(int j=0;j<k;++j)for(int c=0;c<6;++c) require(history.derivatives[j][c]==y[j][c],"History differs from NAO kernel");
            elapsed+=h;require(elapsed==reference_time,"Time differs from original kernel");
        }
        // 原始权重与多项式积分矩条件，独立于轨道力模型。
        for(int l=0;l<k;++l) {
            require(history.predictor[l]==model.cache_AdamsCoef.bex[l],"Predictor coefficients changed");
            require(history.corrector[l]==model.cache_AdamsCoef.bim[l],"Corrector coefficients changed");
        }
        for(int p=0;p<k;++p) {
            long double pred=0,corr=0;
            for(int l=0;l<k;++l) {
                pred+=history.predictor[l]*std::pow(static_cast<long double>(-l),p);
                corr+=history.corrector[l]*std::pow(static_cast<long double>(1-l),p);
            }
            require(std::abs(pred-1.L/(p+1))<1e-6L,"Predictor polynomial moment failed");
            require(std::abs(corr-1.L/(p+1))<1e-6L,"Corrector polynomial moment failed");
        }
    }
    std::cout<<"Original NAO kernel matches: all k=2..9, PE/PECE, both directions, 20 normal steps.\n";
}

void check_failure_reset(OrbitPropagationJob job,const State& initial) {
    const double epoch=UTC2TT(job.opn_setting.t0_MJD_UTC);
    job.opn_setting.adams={3,true};
    State x=initial;AdamsHistory history;double elapsed=0,accepted,error;
    bastion::Overwatch setup(__FILE__,"setup_history",__LINE__);
    for(int i=0;i<4;++i){Adams_NSG2(setup,job,epoch,elapsed,x.data(),.01,history,accepted,error);elapsed+=accepted;}
    require(setup.Code()==ERR_NOERROR && history.filled==3,"Cannot set up history");
    auto fails=[&](OrbitPropagationJob bad,double t0,double t,double h,State input,const char* message) {
        const State original=input;AdamsHistory cache=history;
        bastion::Overwatch trace(__FILE__,"expected_failure",__LINE__);
        Adams_NSG2(trace,bad,t0,t,input.data(),h,cache,accepted,error);
        require(trace.Code()!=ERR_NOERROR && accepted==0 && std::isnan(error),message);
        for(int c=0;c<6;++c) require(input[c]==original[c] || (std::isnan(input[c]) && std::isnan(original[c])),"Failed step changed state");
        require(equalHistory(cache,history),"Failed step changed history");
        std::ostringstream stack,events;trace.ReportTrace(stack);trace.ReportEvents(events);
        require(stack.str().find("Adams_NSG2")!=std::string::npos && events.str().find("[ERROR]")!=std::string::npos,"Failure lost tracer context");
    };
    for(int k:{0,1,10}){auto bad=job;bad.opn_setting.adams.steps=k;fails(bad,epoch,elapsed,.01,x,"Invalid k accepted");}
    for(double h:{0.,1e-20,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        fails(job,epoch,elapsed,h,x,"Invalid step accepted");
    for(double bad:{-1.,0.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid=job;invalid.opn_setting.rkf.maxTE_norm=bad;
        fails(invalid,epoch,elapsed,.01,x,"Invalid startup tolerance accepted");
    }
    for(double bad:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        fails(job,bad,elapsed,.01,x,"Nonfinite epoch accepted");
        fails(job,epoch,bad,.01,x,"Nonfinite elapsed time accepted");
        State invalid=x;invalid[0]=bad;fails(job,epoch,elapsed,.01,invalid,"Nonfinite state accepted");
    }
    {
        AdamsHistory cache=history; bastion::Overwatch trace(__FILE__,"null_state",__LINE__);
        Adams_NSG2(trace,job,epoch,elapsed,nullptr,.01,cache,accepted,error);
        require(trace.Code()!=ERR_NOERROR && accepted==0 && equalHistory(cache,history),"Null state not rejected safely");
    }
    auto missing=job;missing.opn_EOP=EOP{};
    fails(missing,epoch,elapsed,.01,x,"Missing EOP accepted");
    // EOP 在 60002 日结束。通过旧历元的起点评估，但末端评估失败，仍须完整回滚。
    fails(job,60001.99,0,10.,x,"Endpoint EOP failure accepted");
    auto strict=job;strict.opn_setting.rkf.maxTE_norm=1e-30;strict.opn_setting.rkf.iterCount_max=1;
    fails(strict,epoch,elapsed,.5,x,"Startup error limit ignored");
    strict.opn_setting.rkf.fixed_step=true;
    fails(strict,epoch,elapsed,.5,x,"Fixed RKF bypassed Adams startup tolerance");
    strict.opn_setting.rkf.fixed_step=false;
    strict.opn_setting.rkf.maxTE_norm=1e-12;strict.opn_setting.rkf.iterCount_max=30;
    fails(strict,epoch,elapsed,2.,x,"Shortened startup step accepted into fixed history");
    // 步长、方向、步数、模式、状态与时间不连续都必须重新起步。
    for(int mode=0;mode<6;++mode) {
        auto changed=job;State state=x;AdamsHistory cache=history;double h=.01,t=elapsed;
        if(mode==0)h=.005;if(mode==1)h=-.01;
        if(mode==2)changed.opn_setting.adams.steps=4;
        if(mode==3)changed.opn_setting.adams.use_correction=false;
        if(mode==4)state[0]+=1e-6;if(mode==5)t+=.01;
        bastion::Overwatch trace(__FILE__,"reset",__LINE__);
        Adams_NSG2(trace,changed,epoch,t,state.data(),h,cache,accepted,error);
        require(trace.Code()==ERR_NOERROR && cache.filled==1 && std::isfinite(error),"Discontinuous call reused stale history");
    }
    history.Reset();require(history.filled==0 && history.steps==0,"Explicit Reset failed");
}

void check_analytic_and_endpoints(OrbitPropagationJob job) {
    job.opn_model=OrbitPropagationNModel{};
    double max_pos=0,max_vel=0;
    for(int k:{2,5,9})for(bool correction:{false,true})for(double duration:{7.13,-7.13,.03,-.03,.003,-.003,0.}) {
        job.opn_setting.adams={k,correction};job.opn_setting.duration_norm=duration;job.opn_setting.h0_norm=.01;
        std::vector<double> epochs;std::vector<NonSingular2> states;
        bastion::Overwatch trace(__FILE__,"analytic",__LINE__);
        OrbitPropagation_Adams_NSG2(trace,job,epochs,states);
        require(trace.Code()==ERR_NOERROR && !epochs.empty() && epochs.size()==states.size(),"Analytic propagation failed");
        double elapsed=0;
        for(size_t i=0;i<states.size();++i) {
            if(i) elapsed=std::copysign(std::min(i*.01,std::abs(duration)),duration);
            State cart{};states[i].ToCartesian(trace,cart.data());
            State expected=analytic(elapsed);
            max_pos=std::max(max_pos,separation(cart,expected,0)*constants::Scale::UNIT_LENGTH_M);
            max_vel=std::max(max_vel,separation(cart,expected,3)*constants::Scale::UNIT_LENGTH_M/constants::Scale::UNIT_TIME_S);
            if(i)require((epochs[i]-epochs[i-1])*duration>0,"Epoch direction wrong");
        }
        double expected=duration==0?job.opn_setting.t0_MJD_UTC:TT2UTC(UTC2TT(job.opn_setting.t0_MJD_UTC)+duration*constants::Conversion::NORMALIZED_TIME_TO_DAY);
        require(std::abs(epochs.back()-expected)*86400<1e-6,"Partial final step endpoint wrong");
        job.opn_setting.h0_norm=-.01;
        std::vector<double> e2;std::vector<NonSingular2>s2;bastion::Overwatch t2(__FILE__,"step_sign",__LINE__);
        OrbitPropagation_Adams_NSG2(t2,job,e2,s2);
        std::array<double, 6> exported_values_1; states.back().ToArray(exported_values_1.data());
        std::array<double, 6> exported_values_2; s2.back().ToArray(exported_values_2.data());
        require(e2==epochs && exported_values_2==exported_values_1,"Duration did not control direction");
    }
    require(max_pos<1e-4 && max_vel<1e-7,"Two-body analytic accuracy exceeded limit");
    std::cout<<"METRIC two_body_max_position_m="<<max_pos<<" velocity_mps="<<max_vel<<'\n';
    for(double invalid:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        auto bad=job;bad.opn_setting.t0_MJD_UTC=invalid;
        std::vector<double> e;std::vector<NonSingular2>s;bastion::Overwatch t(__FILE__,"invalid_utc",__LINE__);
        OrbitPropagation_Adams_NSG2(t,bad,e,s);
        require(t.Code()!=ERR_NOERROR && e.empty() && s.empty(),"Nonfinite UTC reached time conversion/output");
    }
    // 后期力模型失败不应清空已接受记录，也不能追加坏状态。
    job.opn_setting.t0_MJD_UTC=60001.998;job.opn_setting.duration_norm=1.;job.opn_setting.h0_norm=.01;
    std::vector<double> e;std::vector<NonSingular2>s;bastion::Overwatch trace(__FILE__,"partial_failure",__LINE__);
    OrbitPropagation_Adams_NSG2(trace,job,e,s);
    require(trace.Code()!=ERR_NOERROR && s.size()>1 && s.size()==e.size() && e.back()<60002.,"Partial ephemeris failure not protected");
}
}
int main() try {
    namespace fs=std::filesystem;
    fs::path root=fs::temp_directory_path()/("cdc-adams-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(root);
    struct Cleanup{fs::path p;~Cleanup(){std::error_code ec;fs::remove_all(p,ec);}} cleanup{root};
    fs::path file=root/"eop.txt";std::ofstream out(file);
    out<<"BEGIN OBSERVED\n2023 2 24 59999 .1 .2 .3 0 0 0 0 0 37\n2023 2 25 60000 .1 .2 .3 0 0 0 0 0 37\n2023 2 26 60001 .1 .2 .3 0 0 0 0 0 37\n2023 2 27 60002 .1 .2 .3 0 0 0 0 0 37\nEND OBSERVED\n";out.close();
    OrbitPropagationJob job{};std::string filename=file.string();bastion::Overwatch setup(__FILE__,"setup",__LINE__);
    job.opn_EOP.Read_EOP(setup,filename);require(setup.Code()==ERR_NOERROR,"EOP setup failed");
    job.opn_setting.t0_MJD_UTC=60000.25;job.opn_setting.display_progress=0;job.opn_setting.rkf.iterCount_max=30;
    State cart=analytic(0);job.initCart_norm.Set(cart[0],cart[1],cart[2],cart[3],cart[4],cart[5],true);
    job.initCart_norm.ToKeplerian(setup,job.initKep_norm);
    NonSingular2 nsg;job.initCart_norm.ToNonSingular2(setup,nsg);
    require(setup.Code()==ERR_NOERROR,"Orbit setup failed");
    std::cout<<std::setprecision(17);
    std::array<double, 6> exported_values_3; nsg.ToArray(exported_values_3.data());
    check_reference(job,exported_values_3);
    std::array<double, 6> exported_values_4; nsg.ToArray(exported_values_4.data());
    check_failure_reset(job,exported_values_4);
    check_analytic_and_endpoints(job);
    std::cout<<"Adams checks passed: "<<checks<<'\n';return 0;
} catch(const std::exception& e){std::cerr<<"Adams test failed: "<<e.what()<<'\n';return 1;}
