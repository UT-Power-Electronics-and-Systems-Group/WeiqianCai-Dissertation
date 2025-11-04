clear;clc;
% System Parameters
w0 = 2*pi*50;
Lg = 7e-6+3.5e-6+100*31.69e-6; %*100;% Lg, Line, Lxfmr
vth = 1.3e3;
vdcref = 1.5e3; 
mp = pi/125000; % 0.5Hz, 125kVA
 
Pref = 50e3; 
wc = 2*pi*100;   
V0 = 600*sqrt(2)/sqrt(3); % L-2-N PEAK   
Xg = w0*Lg;
Cdc = 4*1130e-6; 
Pmpp = 120e3;

% Time span
tspan = [0 2]; 

%============== Kpmax ==============
kp = pi/(vdcref-vth);  

DELTA = asin(2*Xg*Pmpp/(3*V0^2));
DELTA_approx = 2*Xg*Pmpp/(3*V0^2);
VDC = mp/kp*(Pmpp-Pref) + vdcref
PF = Pmpp;
Gpd = (3*V0^2*cos(DELTA))/(2*Xg);
k_m = Gpd*mp*mp*Cdc/4;
k_b = k_m * vdcref;
k_c = k_m * mp * (Pmpp-Pref);
k_d = sqrt(k_b*k_b+4*k_c);
kp_min = 1.0 * mp*Pmpp/(0.5*vdcref)
kp_max = 0.5*(k_b+k_d)
kesi = 0.5*Gpd*mp/sqrt(Gpd*kp/(Cdc*VDC))

% Parameter vector
Params = [mp, kp, Pref, vdcref, wc, V0, Xg, Cdc, Pmpp];

% Initial States for simplified model
P0 = 140e3;
vdc0 = vth;
x0_simp = [P0, vdc0];  % Simplified system initial state [P, vdc]

% Solve simplified model only
[t_simp, x_simp] = ode45(@(t, x) fun_linsimp(x, Params), tspan, x0_simp);

% Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% Color scheme
green = [0.4660, 0.6740, 0.1880];

% Power plot
subplot(2,1,1);
plot(t_simp, x_simp(:,1), 'Color', green, 'LineWidth', 2);
grid on;
title('Simplified Model: Active Power $P_{\mathrm{f}}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Power (W)', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'TickLabelInterpreter', 'latex');

% DC voltage plot
subplot(2,1,2);
plot(t_simp, x_simp(:,2), 'Color', green, 'LineWidth', 2);
grid on;
title('Simplified Model: DC Voltage $v_{\mathrm{dc}}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$v_{\mathrm{dc}}$ (V)', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'TickLabelInterpreter', 'latex');

% Formatting adjustments
set(gcf, 'Position', [100, 100, 800, 600]);
h = findobj(gcf, 'type', 'axes');
set(h, 'Box', 'on', 'FontSize', 10);

% Simplified 2nd-order linear model
function dx = fun_linsimp(x, param)  
    % Simplified 2nd-order linear model
    mp = param(1);    kp = param(2);   Pref = param(3); 
    vdcref = param(4); wc = param(5);    V0 = param(6); 
    Xg = param(7);   Cdc = param(8);   Pmpp = param(9); 
    
    DELTA = asin(2*Xg*Pmpp/(3*V0^2));
    PF = Pmpp;
    VDC = mp/kp*(Pmpp-Pref) + vdcref;
    Gpd = (3*V0^2*cos(DELTA))/(2*Xg);

    P = x(1); 
    vdc = x(2); 
    
    dx = zeros(2,1);
    dx(1) = -Gpd*mp*(P-PF) + Gpd*kp*(vdc-VDC);
    dx(2) = -(P-Pmpp)/(Cdc*VDC);
end