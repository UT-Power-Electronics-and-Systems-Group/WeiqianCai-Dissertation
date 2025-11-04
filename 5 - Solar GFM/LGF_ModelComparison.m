clear;clc;
% System Parameters

w0 = 2*pi*50;

Lg = 7e-6+3.5e-6+100*31.69e-6; % Lg, Line, Lxfmr

vth = 1.3e3;
vdcref = 1.5e3; 
mp = pi/125000; % 0.5Hz, 125kVA

kp = pi/(vdcref-vth)
Pref = 50e3; 
wc = 2*pi*30;%100;    
V0 = 600*sqrt(2)/sqrt(3); % L-2-N PEAK   
Xg = w0*Lg;
Cdc = 4*1130e-6; 
Pmpp = 120e3;
np= 3*V0*V0/(2*Xg)*mp;
zeta= sqrt(wc/np)/2
wnp=sqrt(np*wc)


% Time span
tspan = [0 2]; 

%============== Kpmax ==============
DELTA = asin(2*Xg*Pmpp/(3*V0^2));
VDC = mp/kp*(Pmpp-Pref) + vdcref
PF = Pmpp;
Gpd = (3*V0^2*cos(DELTA))/(2*Xg)
k_m = Gpd*mp*mp*Cdc/4;
k_b = k_m * vdcref;
k_c = k_m * mp * (Pmpp-Pref);
k_d = sqrt(k_b*k_b+4*k_c);
kp_min = 1.0 * mp*Pmpp/(0.5*vdcref)
kp_max = 0.5*(k_b+k_d)
kesi = 0.5*Gpd*mp/sqrt(Gpd*kp/(Cdc*VDC))

% Parameter vector
Params = [mp, kp, Pref, vdcref, wc, V0, Xg, Cdc, Pmpp];
Pdelta_ratio = 3*V0*V0/(2*Xg);

% Initial States
P0 = 140e3;
Pf0 = P0;
Delta0 = asin((2*P0*Xg)/(3*V0*V0)); 
vdc0 = vth;
x0 = [Delta0, Pf0, vdc0];          % Original system initial state
x0_simp = [Pf0, vdc0];            % Simplified system initial state


% Solve ODE systems
[t_lin, x_lin] = ode45(@(t, x) fun_linear(x, Params), tspan, x0);
[t_nlin, x_nlin] = ode45(@(t, x) fun_nonlinear(x, Params), tspan, x0);
[t_simp, x_simp] = ode45(@(t, x) fun_linsimp(x, Params), tspan, x0_simp);

% Calculate power signals
P_lin = x_lin(:,1) * Pdelta_ratio;
P_nlin = x_nlin(:,1) * Pdelta_ratio;

% Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 800]);

% Color scheme
blue = [0, 0.4470, 0.7410];
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];

% Delta comparison
subplot(3,1,1);
plot(t_lin, x_lin(:,1), 'Color', blue, 'LineWidth', 1.5);
hold on;
plot(t_nlin, x_nlin(:,1), 'Color', red, 'LineWidth', 1.5, 'LineStyle','--');
hold off;
grid on;
title('State Comparison: $\delta(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\delta$ (rad)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'Linear','Nonlinear'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'TickLabelInterpreter', 'latex');

% Power comparison
subplot(3,1,2);
p1 = plot(t_lin, x_lin(:,2), 'Color', blue, 'LineWidth', 2);
hold on;
p2 = plot(t_nlin, x_nlin(:,2), 'Color', red, 'LineWidth', 2, 'LineStyle','--');
p3 = plot(t_lin, P_lin, 'Color', blue, 'LineWidth', 1.5, 'LineStyle',':');
p4 = plot(t_nlin, P_nlin, 'Color', red, 'LineWidth', 1.5, 'LineStyle','-.');
p5 = plot(t_simp, x_simp(:,1), 'Color', green, 'LineWidth', 2, 'LineStyle','-');
grid on;
title('Power Comparison: $P_{\mathrm{f}}(t)$ vs $P(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Power (W)', 'Interpreter', 'latex', 'FontSize', 10);
legend([p1,p2,p3,p4,p5],...
    {'$P_{\mathrm{f,linear}}$','$P_{\mathrm{f,nonlinear}}$',...
    '$P_{\mathrm{linear}}$','$P_{\mathrm{nonlinear}}$',...
    '$P_{\mathrm{simplified}}$'},...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'TickLabelInterpreter', 'latex');

% DC voltage comparison
subplot(3,1,3);
plot(t_lin, x_lin(:,3), 'Color', blue, 'LineWidth', 1.5);
hold on;
plot(t_nlin, x_nlin(:,3), 'Color', red, 'LineWidth', 1.5, 'LineStyle','--');
plot(t_simp, x_simp(:,2), 'Color', green, 'LineWidth', 1.5, 'LineStyle','-');
hold off;
grid on;
title('State Comparison: $v_{\mathrm{dc}}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$v_{\mathrm{dc}}$ (V)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'Linear','Nonlinear','Simplified'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'TickLabelInterpreter', 'latex');

% Formatting adjustments
set(gcf, 'Position', [100, 100, 800, 800]);
h = findobj(gcf, 'type', 'axes');
set(h, 'Box', 'on', 'FontSize', 10);

%% Models: original nonlinear model, 3-order linear model, 2nd-order linear model
% ----- System Models -----
function dx = fun_nonlinear(x, param)  
    % Original nonlinear model
    mp = param(1);    kp = param(2);   Pref = param(3); 
    vdcref = param(4); wc = param(5);    V0 = param(6); 
    Xg = param(7);   Cdc = param(8);   Pmpp = param(9); 
    
    delta = x(1); 
    Pf = x(2); 
    vdc = x(3); 
    
    dx = zeros(3,1);
    dx(1) = mp*(Pref-Pf) + kp*(vdc-vdcref);
    dx(2) = -wc*Pf + (3*wc*V0^2*sin(delta))/(2*Xg);
    dx(3) = (Pmpp - (3*V0^2*sin(delta))/(2*Xg))/(Cdc*vdc);
end

function dx = fun_linear(x, param)  
    % Original linearized model
    mp = param(1);    kp = param(2);   Pref = param(3); 
    vdcref = param(4); wc = param(5);    V0 = param(6); 
    Xg = param(7);   Cdc = param(8);   Pmpp = param(9); 
    
    DELTA = asin(2*Xg*Pmpp/(3*V0^2));
    PF = Pmpp;
    VDC = mp/kp*(Pmpp-Pref) + vdcref;

    delta = x(1); 
    Pf = x(2); 
    vdc = x(3); 
    
    dx = zeros(3,1);
    dx(1) = -mp*(Pf-PF) + kp*(vdc-VDC);
    dx(2) = (3*wc*V0^2*cos(DELTA))/(2*Xg)*(delta-DELTA) - wc*(Pf-PF);
    dx(3) = -(3*V0^2*cos(DELTA))/(2*Xg*Cdc*VDC)*(delta-DELTA);
end

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
    dx(1) = -Gpd*mp      * (P-PF) + Gpd*kp * (vdc-VDC);
    dx(2) = -1/(Cdc*VDC) * (P-PF);
end
