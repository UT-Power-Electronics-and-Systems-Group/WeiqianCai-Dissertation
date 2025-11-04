clear;clc;
%% PV Boost Params & Curve Generation
Vdc = 1500;
fsw = 50e3;
wsw = 2*pi*fsw;
% global Voc Isc Vmpp Impp
Voc = 1400;%1400;
Isc = 40;%160;
Vmpp = 1000;%1000;
Impp = 30;%120;
Rmpp = Vmpp/Impp; 
Pmpp = Vmpp*Impp;
Vin = Vmpp;
Dn = Vin/Vdc;
D  = 1-Dn;
Lmin = D*Rmpp/(0.6*fsw)
Cmin = D*(Pmpp/Vdc)/(5*fsw)

L = 400e-6;%120e-6;
rc = 1.4e-3;%0.1;
C = 1130e-6;%120e-6;
Ci= 14e-6;%56e-6;

Rload = Vdc^2/30e3;%120e3;

% Simplified PV Curve
% [Ref] Bellini, Armando, et al. "Simplified model of a photovoltaic module." 2009 Applied Electronics. IEEE, 2009.
C2_PV = (Vmpp/Voc-1) / log(1-Impp/Isc);
C1_PV = (1-Impp/Isc) * exp( -Vmpp/(C2_PV*Voc) );

% Temperature Coeff
% Alpha_PV = 1.46e-3*Isc/5.3;
% Beta_PV  = -158e-3*Voc/44.6;

% Irradiance Coeff - G_PV, s.t.
% Isc(G_PV) = Iscs  * G_PV/Gs_PV; % Iscs @ Gs_PV=1kW/m^2, Ts=25 degC.
% Impp(G_PV)= Impps * G_PV/Gs_PV;

Vpv = 0:1400;
Ipv = Isc .* ( 1 - C1_PV .* ( exp( Vpv./(C2_PV*Voc) ) - 1 ) );
Ppv = Vpv .* Ipv;

plot(Vpv, Ipv);
% hold on
% plot(Vpv, Ppv);
%% Boost Modeling
wzi = 2/(Rload*C);
wzv = Dn^2*Rload/(L);
wesr = 1/(rc*C);
wo  = Dn/sqrt(L*C);
To  = 1/wo;
Q   = Dn*Rload*sqrt(C/L);
fprintf('fsw=%f, fesr=%f, fzi=%f, fzv=%f, fo=%f\n', ...
         wsw/(2*pi),wesr/(2*pi),wzi/(2*pi),wzv/(2*pi),wo/(2*pi));
% freq_pz_kHz=[wsw,wesr,wzv,wo]/(2*pi*1000)

s= tf('s');
opts = bodeoptions('cstprefs');
% SINCE wesr>> curr. loop BW, even >> fs, so it can be neglected
% Gvi = (Dn*Rload/2)* (1-s/wzv)*(1+s/wesr)/(1+s/wzi);
Gvi = (Dn*Rload/2)* (1-s/wzv)/(1+s/wzi);
Gid = (Vdc/L) * (s+wzi) / ( s^2 + (wo/Q)*s + wo^2 );
% control-to-output w/o cap. ESR
Gvd = (Vdc/Dn) * (1-s/wzv)/( (s/wo)^2 + s/(Q*wo) + 1 );
% control-to-output with cap. ESR
Gvdesr = (Vdc/Dn) * (1-s/wzv)*(1+s/wesr)/( (s/wo)^2 + s/(Q*wo) + 1 );
% 1.5 Tctrl delay - not considered in voltage. control
Gdl  = 1/(1+1.5*s*(1/fsw));

%% Double Loop - Current Compensator
%When kp=1, 20lg|Gid*Gdl|=-81.1dB --> kp=11350, BW=5kHz
Gcc = 481*(1/s)*(1+s/wo)^2/(s/wzi+1); 
Gp_i = Gid*Gdl;
Gol_i= Gid*Gdl*Gcc;
Gcl_i = Gol_i/(1+Gol_i);

% Curr Loop Bode
% bodeplot(Gid)
% hold on
bodeplot(Gp_i)
hold on 
bodeplot(Gol_i)
hold on
bodeplot(Gcl_i)

%% Double Loop - Voltage Compensator
% wcv=2*pi*500;
Gcv = 0.06659 + 2.885/s;
Gvi_p=Gvi*Gcl_i;
Gol_v=Gcv*Gvi;
Gol_v2=Gcv*Gvi_p;
Gcl_v = Gol_v2/(Gol_v2+1);

% bodeplot(Gvi)
% hold on
bodeplot(Gvi_p)
hold on 
bodeplot(Gol_v2)
% hold on
% bodeplot(Gcl_v)

%% [Inactive] Design V2 with small LC & single voltage compensator
% BW_max = [wzv/(5*2*pi), 0.1*fsw] % fcv<fzv/5, fcv<fsw/10
% BW_min = wo/pi % fcv > 2*wo
% Gc_leadlag = (1/s)*(1+s/wo)^2/( (1+s/wesr) * (1+s/wzv) )*3.7584;
% Gc_leadlag = (1/s)*(1+s/wo)^2/( (1+s*2/wsw) * (1+s/wzv) );
% h=bodeplot(Gvdesr*Gc_leadlag);
% h=bodeplot(Gvdesr);
% hold on
% bodeplot(Gvdesr*Gc_leadlag)

%% [Inactive] Boost Control Design considering PV model
% vpv = 750:5:1200;
% vpv = [800,900,1000,1100];
% len=length(vpv);
% % Linear I/V at CVR
% % ipv = Impp - Impp/(Voc-Vmpp)*(vpv-Vmpp);
% % rpv = (Voc-Vmpp)/Impp.*ones(size(vpv));
% % Ppv = vpv.*ipv;
% % Quqadratic P/V at CVR
% K1 = Impp*Vmpp/((Voc-Vmpp)^2);
% K2 = Impp*Vmpp;
% Ppv = K2 - K1.*(vpv-Vmpp).^2;
% ipv = Ppv./vpv;
% rpv = 1./(K1-(K1*Vmpp*Vmpp-K2)./(vpv.^2));
% 
% Rpv = vpv./ipv;
% pv_list = [vpv;ipv;Ppv;Rpv;rpv];
% Dn = vpv/Vdc; % Dn = 1-D

% wz2_0 = Rpv./L - 1./(rpv.*Ci)
% wz2_1 = Rpv./L
% s=tf('s');
% for i=1:len
% Gvd_0(i) = -ipv(i)/C...
% *( s^2 - s*( Rpv(i)/L -1/(rpv(i)*Ci) ) + (1-Rpv(i)/rpv(i))/(L*C))...
% /( s^3 + s^2/(rpv(i)*Ci) + s*(C+Dn(i)^2*Ci)/(L*C*Ci) + Dn(i)^2/(L*C*Ci*rpv(i)) )
% % bodeplot(Gvd_0(i))
% Gvd_1(i) = -ipv(i)/C...
% * (s-(1-Rpv(i)/rpv(i))/(Rpv(i)*Ci)) * (s-Rpv(i)/L)...
% / ((s+Dn(i)^2/(rpv(i)*C)) * (s^2+1/(rpv(i)*Ci)*s+1/(L*Ci)) )
% hold on
% end

% bodeplot(Gvd_0(3))
% hold on
% bodeplot(Gvd_1(3))
% hold off

% plot(vpv,ipv/Impp);
% hold on
% plot(vpv,Rpv/Rmpp);
% hold on
% plot(vpv,rpv/Rmpp);

