#include <array>
#include "orbit_expansion.hpp"
#include "analyt_averaging.hpp"
#include "model_ephemeris_analytic.hpp"
#include "orbit_statedef.hpp"
#include "math_basic.hpp"
#include "utils_constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
void near(double x,double y,double tolerance,const char* message) {
    require(std::isfinite(x)&&std::abs(x-y)<=tolerance,message);
}
orbit::NonSingular2 state(bastion::Overwatch& trace,double inclination=0) {
    orbit::Keplerian kep(7000000./ae,.001,inclination,0.,.4,.6,true);
    orbit::NonSingular2 nsg;
    kep.ToNonSingular2(trace,nsg);
    return nsg;
}
}

int main() try {
    // 独立低阶封闭式同时核对倾角函数的符号、值和两阶导数。
    for (double i:{0.,1e-9,.1,.7,1.5,3.0,pi}) {
        double F,dF,ddF;
        FlmpDerivatives(2,0,1,i,F,dF,ddF);
        near(F,.75*pow(sin(i),2)-.5,2e-14,"F201");
        near(dF,1.5*sin(i)*cos(i),2e-14,"F201 derivative");
        near(ddF,1.5*cos(2*i),2e-14,"F201 second derivative");
        FlmpDerivatives(2,1,0,i,F,dF,ddF);
        near(F,.75*sin(i)*(1+cos(i)),2e-14,"F210");
        near(dF,.75*(cos(i)+cos(2*i)),2e-14,"F210 derivative including zero inclination");
        near(ddF,-.75*(sin(i)+2*sin(2*i)),2e-14,"F210 second derivative");
        FlmpDerivatives(2,2,0,i,F,dF,ddF);
        near(F,.75*pow(1+cos(i),2),2e-14,"F220");
    }
    // 高阶导数用五点数值差分独立核对；尺度取系数和量级，避免过零处相对误差失真。
    int derivative_checks=0;
    for (int l=2;l<=10;++l) for (int m=0;m<=l;++m) for (int p=0;p<=l;++p) {
        for (double i:{.1,.7,1.5,2.8}) {
            double F,dF,ddF;
            FlmpDerivatives(l,m,p,i,F,dF,ddF);
            const double h=1e-4;
            const double fm2=Flmp(l,m,p,i-2*h),fm1=Flmp(l,m,p,i-h);
            const double fp1=Flmp(l,m,p,i+h),fp2=Flmp(l,m,p,i+2*h);
            const double d=(-fp2+8*fp1-8*fm1+fm2)/(12*h);
            const double dd=(-fp2+16*fp1-30*F+16*fm1-fm2)/(12*h*h);
            const double scale=std::max({1.,std::abs(F),std::abs(dF),std::abs(ddF)});
            near(dF,d,2e-7*scale,"Inclination first derivative finite difference");
            near(ddF,dd,3e-5*scale,"Inclination second derivative finite difference");
            ++derivative_checks;
        }
        for (double i:{0.,1e-10,pi-1e-10,pi}) {
            double F,dF,ddF; FlmpDerivatives(l,m,p,i,F,dF,ddF);
            require(std::isfinite(F)&&std::isfinite(dF)&&std::isfinite(ddF),"Nonfinite endpoint inclination function");
        }
    }
    double X,dX;
    HansenSeries(-3,0,0,.1,4,X,dX);
    near(X,1+1.5*.01+1.875*.0001,1e-15,"Hansen zero harmonic");
    near(dX,.3,1e-15,"Hansen truncated derivative");
    HansenSeries(-3,0,1,.1,4,X,dX);
    near(X,.15+27./16*.001,1e-15,"Hansen first harmonic");
    HansenSeries(-3,0,0,0,4,X,dX);
    near(X,1,0,"Circular Hansen value"); near(dX,0,0,"Circular Hansen derivative");

    bastion::Overwatch trace(__FILE__,__func__,__LINE__);
    // 日月分别核对新旧调用：可选辅助速度不能改变原三分量位置，包括最后一位。
    for (bool moon : {false,true}) {
        for (double epoch : {51544.5,60149.}) {
            double position_only[3],position[3],auxiliary_velocity[3];
            if (moon) {
                Cart_MeanMoon(epoch,position_only);
                Cart_MeanMoon(epoch,position,auxiliary_velocity);
            } else {
                Cart_MeanSun(epoch,position_only);
                Cart_MeanSun(epoch,position,auxiliary_velocity);
            }
            // 不硬编码另一套历表答案，只要求同一次原公式的位置逐分量完全相同。
            for (int j=0;j<3;++j) {
                require(position_only[j]==position[j],
                    moon ? "Optional Moon velocity changed position" : "Optional Sun velocity changed position");
            }

            // 辅助速度只用于 mu=1 下恢复方向；转根数后重构位置，不把它当作真实日月速度。
            orbit::Cartesian cart(position[0],position[1],position[2],
                auxiliary_velocity[0],auxiliary_velocity[1],auxiliary_velocity[2],true);
            orbit::Keplerian elements;
            cart.ToKeplerian(trace,elements);
            require(trace.Code()==bastion::ERR_NOERROR,
                moon ? "Moon auxiliary state to elements failed" : "Sun auxiliary state to elements failed");

            // 采用当前位置量级的浮点容差；源位置仍来自原函数，不另存参考数值。
            double rebuilt[6];
            elements.ToCartesian(trace,rebuilt);
            require(trace.Code()==bastion::ERR_NOERROR,
                moon ? "Moon auxiliary elements reconstruction failed" : "Sun auxiliary elements reconstruction failed");
            const double position_scale=std::max({1.,std::abs(position[0]),std::abs(position[1]),std::abs(position[2])});
            for (int j=0;j<3;++j) {
                near(rebuilt[j],position[j],2e-12*position_scale,
                    moon ? "Moon auxiliary state position reconstruction" : "Sun auxiliary state position reconstruction");
            }
        }
    }

    orbit::AveragingOptions opt;
    const auto input=state(trace,.3);
    auto output=input;
    opt.use_j2=opt.use_j22=opt.use_frame=opt.use_jl=opt.use_jlm=false;
    opt.use_lunar=opt.use_solar=opt.use_srp=false;
    orbit::OsculatingToMean(trace,60149.,input,opt,output);
    require(trace.Code()==bastion::ERR_NOERROR,"Identity conversion failed");
    std::array<double, 6> exported_values_1; output.ToArray(exported_values_1.data());
    std::array<double, 6> exported_values_2; input.ToArray(exported_values_2.data());
    const auto a=exported_values_2,b=exported_values_1;
    for (int j=0;j<6;++j) near(a[j],b[j],1e-15,"All-off NSG2 identity");

    // JL 单独启用也必须更新 F/G；参考实现依赖 J22/JLM 开关的缓存条件不适用于本库。
    opt.use_jl=true;
    orbit::OsculatingToMean(trace,60149.,input,opt,output);
    std::array<double, 6> exported_values_3; output.ToArray(exported_values_3.data());
    require(trace.Code()==bastion::ERR_NOERROR && exported_values_3!=a,"Standalone JL produced no correction");
    opt=orbit::AveragingOptions{};
    auto zero=state(trace),near_zero=state(trace,1e-10);
    orbit::NonSingular2 zero_mean,near_mean;
    orbit::OsculatingToMean(trace,60149.,zero,opt,zero_mean);
    orbit::OsculatingToMean(trace,60149.,near_zero,opt,near_mean);
    require(trace.Code()==bastion::ERR_NOERROR,"Zero inclination conversion failed");
    std::array<double, 6> exported_values_4; near_mean.ToArray(exported_values_4.data());
    std::array<double, 6> exported_values_5; zero_mean.ToArray(exported_values_5.data());
    auto z0=exported_values_5,z1=exported_values_4;
    for (int j=0;j<6;++j) near(z0[j],z1[j],1e-8,"O2M discontinuity at zero inclination");

    // 公共接口失败时不改输出，包括超出解析范围和非法选项。
    for (int kind=0;kind<3;++kind) {
        bastion::Overwatch failed(__FILE__,__func__,__LINE__);
        auto invalid=input; auto result=input; auto badopt=opt;
        if (kind==0) badopt.hansen_order=5;
        if (kind==1) invalid.Set(failed,100000000./ae,.01,0,0,0,0,true);
        const double epoch=kind==2 ? std::numeric_limits<double>::quiet_NaN() : 60149.;
        orbit::OsculatingToMean(failed,epoch,invalid,badopt,result);
        require(failed.Code()!=bastion::ERR_NOERROR,"Invalid averaging call succeeded");
        std::array<double, 6> exported_values_6; result.ToArray(exported_values_6.data());
        require(exported_values_6==a,"Failed averaging call modified output");
    }
    std::cout<<"Orbit averaging checks passed; "<<derivative_checks<<" derivative checks\n";
    return 0;
} catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
