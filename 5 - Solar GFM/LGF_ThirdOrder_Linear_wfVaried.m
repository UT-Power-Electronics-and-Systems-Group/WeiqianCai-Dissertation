%%
clear;clc;
% System Parameters
w0 = 2*pi*50;
Lg = 7e-6+3.5e-6+100*31.69e-6; % Lg, Line, Lxfmr
vth = 1.3e3;
vdcref = 1.5e3; 
mp = pi/125000; % 0.5Hz, 125kVA
kp = pi/(vdcref-vth)
Pref = 50e3;   
V0 = 600*sqrt(2)/sqrt(3); % L-2-N PEAK   
Xg = w0*Lg;
Cdc = 4*1130e-6; 
Pmpp = 120e3;


% Steady States & Related Vars
DELTA = asin(2*Xg*Pmpp/(3*V0^2));
VDC = mp/kp*(Pmpp-Pref) + vdcref
PF = Pmpp;
Gpd = (3*V0^2*cos(DELTA))/(2*Xg);
wc_min = 0.5*Gpd*mp

% Initial States
P0 = 140e3;
Pf0 = P0;
Delta0 = asin((2*P0*Xg)/(3*V0*V0)); 
vdc0 = vth;
x0 = [Delta0, Pf0, vdc0];          % Original system initial state

% Time span
tspan = [0 2]; 

% wc sweep
wc_values = wc_min * [2, 5, 10];
num_wc = length(wc_values);

% Simulate again with new wc_values
t_cell = cell(1, num_wc);
x_cell = cell(1, num_wc);
for i = 1:num_wc
    wc = wc_values(i);
    Params = [mp, kp, Pref, vdcref, wc, V0, Xg, Cdc, Pmpp];
    [t_sim, x_sim] = ode45(@(t, x) fun_linear(x, Params), tspan, x0);
    t_cell{i} = t_sim;
    x_cell{i} = x_sim;
end

%% Visualization: 1x3 subplot (2-column width), labeled (a), (b), (c)
figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7.2, 1.7]);

c_pf  = [0, 0.4470, 0.7410];        % blue
c_pac = [0.8500, 0.3250, 0.0980];   % red
nf_subtitles = {'(a) $N_\mathrm{f}=2$', '(b) $N_\mathrm{f}=5$', '(c) $N_\mathrm{f}=10$'};
ax_handles = gobjects(1, num_wc);

% Subplot layout
lefts = [0.08, 0.38, 0.68];  % X positions
width = 0.25;
height = 0.55;               % shorter height
bottom = 0.3;                % higher bottom

for i = 1:num_wc
    ax = axes('Position', [lefts(i), bottom, width, height]);
    ax_handles(i) = ax;

    t = t_cell{i};
    x = x_cell{i};
    Pf = x(:,2) / 1e3;
    Pac = Gpd * x(:,1) / 1e3;

    h1 = plot(t, Pf,  '-',  'Color', c_pf,  'LineWidth', 1.2); hold on;
    h2 = plot(t, Pac, '-', 'Color', c_pac, 'LineWidth', 1.2); hold off;
    grid on;

    xlabel('$t$, [s]', ...
        'Interpreter', 'latex', ...
        'FontSize', 8, ...
        'Color', [0 0 0]);  % Black label
    if i == 1
        ylabel('$P_{\mathrm{f}},\ P_{\mathrm{ac}},$ [kW]', ...
            'Interpreter', 'latex', ...
            'FontSize', 8, ...
            'Color', [0 0 0]);  % Black label
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, ...
        'TickLabelInterpreter', 'latex', ...
        'FontSize', 7, ...
        'XColor', [0 0 0], ...
        'YColor', [0 0 0]);  % Black tick labels
    xlim([0 2]); ylim auto

    % Add subtitle manually as annotation
    pos = get(ax, 'Position');
    center_x = pos(1) + pos(3)/2;
    subtitle_y = pos(2) - 0.08;

    annotation('textbox', ...
        [center_x - 0.08, subtitle_y, 0.16, 0.05], ...
        'String', nf_subtitles{i}, ...
        'Interpreter', 'latex', ...
        'FontSize', 8, ...
        'Color', [0 0 0], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'EdgeColor', 'none');
end

% Legend inside 3rd plot
legend(ax_handles(3), {'$P_{\mathrm{f}}$', '$P_{\mathrm{ac}}$'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 7, ...
    'Location', 'southeast', ...
    'Box', 'on', ...
    'TextColor', [0 0 0], ...
    'EdgeColor', [0 0 0], ...
    'Color', [1 1 1]);

%% Models: 3-order linear model

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