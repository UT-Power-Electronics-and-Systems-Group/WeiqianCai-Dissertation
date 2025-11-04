clear;clc;

Sb=1500;
Vb=120;
Ib=Sb/Vb;
Zb=Vb/Ib;
Yb=Ib/Vb;
w0=2*pi*60;
Lb=Zb/w0;
Cb=1/(w0*Zb);
fsw=100e3;
wsw=2*pi*fsw;
Cf=7e-6;
Lf=300-6;
Rf=0.05;
Lg=3e-3;
Rg=Rf;

kpi=2*pi*4.5e3*Lf;
kii=kpi*100;
kpv=2*pi*2000*Cf;
kiv=kpv*50;

%%
kii_pu=kii/Zb;
kpi_pu=kpi/Zb;
Rf_pu=Rf/Zb;
Lf_pu=Lf/Lb;
wis=kii_pu/(kpi_pu+Rf_pu);
wif=w0/Lf_pu*(kpi_pu+Rf_pu);

Cf_pu=Cf/Cb;
kpv_pu=kpv/Yb;
kiv_pu=kiv/Yb;
wvs=kiv_pu/(kpv_pu);
wvf=w0/Cf_pu*(kpv_pu);
Rc=5.0*2.3127e+04;
Rc_pu=Rc/Zb;
wvs=kiv_pu/(kpv_pu+1/Rc_pu);
wvf=w0/Cf_pu*(kpv_pu+1/Rc_pu);
%% 
wc=2*pi*20;
D=50*Sb/w0;
J=0.2*2*Sb/w0;
mu=2*pi*0.5/(0.9^2*(1-0.9^2)*2*Vb^2);
eta=2*0.5*pi*Vb^2/Sb;
mp=0.5*2*pi/Sb;
mq=0.1*sqrt(2)*Vb/Sb;
mpv=mp*(D+1);
mqv=mq;
mp_pu=0.5*2*pi;
Lg_pu=Lg/Lb;
wgs_droop=mp_pu/Lg_pu;
wgf_droop=wc/2;
wgs_vsm=wgs_droop;
wgf_vsm=(D+1)/(2*J);
eta_pu=eta/Zb;
mu_pu=mu*Vb^2;
wgs_dvoc=eta_pu/Lg_pu
wgf_dvoc=wgf_droop;
[0.0 wis wvs wgs_dvoc wgs_droop wgs_vsm]/(2*pi) %Hz

[wsw wif wvf wgf_dvoc wgf_droop wgf_vsm]/(2*pi*1000) %kHz


