clear; clc;

% System Parameters
w0 = 2*pi*50;
Lg = 7e-6 + 3.5e-6 + 100*31.69e-6; 
vth = 1.3e3;
vdcref = 1.5e3; 
mp = pi/125000; 
Pref = 50e3; 
wc = 2*pi*100;   
V0 = 600*sqrt(2)/sqrt(3); 
Xg = w0*Lg;
Cdc = 4*1130e-6; 
Pmpp = 120e3;

% kp sweep
kp_values = 0.01:0.01:0.03;
num_kp = length(kp_values);

% Initial condition for simplified model
P0 = 140e3;
vdc0 = vth;
x0_simp = [P0, vdc0];
tspan = [0 2];

% Simulate
t_cell = cell(1, num_kp);
x_cell = cell(1, num_kp);
for i = 1:num_kp
    kp = kp_values(i);
    Params = [mp, kp, Pref, vdcref, wc, V0, Xg, Cdc, Pmpp];
    [t_sim, x_sim] = ode45(@(t,x) fun_linsimp(x, Params), tspan, x0_simp);
    t_cell{i} = t_sim;
    x_cell{i} = x_sim;
end

%% Plot
figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 3.5, 2.4]);

cmap = lines(num_kp);
subtitle_labels = {'(a)', '(b)'};

for k = 1:2
    ax = subplot(2,1,k); hold on;
    for i = 1:num_kp
        if k == 1
            ydata = x_cell{i}(:,1) / 1e3;  % P in kW
        else
            ydata = x_cell{i}(:,2);       % Vdc
        end
        plot(t_cell{i}, ydata, 'Color', cmap(i,:), 'LineWidth', 1.1);
    end
    grid on;
    box on;

    % UI element color and font
    set(gca, ...
        'TickLabelInterpreter', 'latex', ...
        'FontSize', 7, ...
        'XColor', [0 0 0], ...
        'YColor', [0 0 0]);  % Set tick label color to black

    if k == 1
        ylabel('$P_{\mathrm{ac}}, P_{\mathrm{f}}$, [kW]', ...
            'Interpreter', 'latex', ...
            'FontSize', 8, ...
            'Color', [0 0 0]);  % Axis label color
    else
        ylabel('$v_{\mathrm{dc}}$, [V]', ...
            'Interpreter', 'latex', ...
            'FontSize', 8, ...
            'Color', [0 0 0]);
        xlabel('$t$, [s]', ...
            'Interpreter', 'latex', ...
            'FontSize', 8, ...
            'Color', [0 0 0]);

        % Add compact legend in 2nd subplot
        lgd = legend(arrayfun(@(kp) sprintf('$\\kappa_\\mathrm{p}=%.2f$', kp), kp_values, ...
            'UniformOutput', false), ...
            'Interpreter', 'latex', ...
            'FontSize', 6, ...
            'Box', 'on', ...
            'Location', 'northeast', ...
            'TextColor', [0 0 0], ...
            'EdgeColor', [0 0 0], ...
            'Color', 'w');
        lgd.ItemTokenSize = [8, 5];  % shorter line
    end

    % Add (a), (b) subtitle at bottom-center
    pos = get(gca, 'Position');
    annotation('textbox', [pos(1), pos(2)-0.06, pos(3), 0.03], ...
        'String', subtitle_labels{k}, ...
        'Interpreter', 'latex', ...
        'FontSize', 7, ...
        'EdgeColor', 'none', ...
        'Color', [0 0 0], ...
        'HorizontalAlignment', 'center');
end


%% Simplified 2nd-order model
function dx = fun_linsimp(x, param)
    mp = param(1); kp = param(2); Pref = param(3);
    vdcref = param(4); wc = param(5); V0 = param(6);
    Xg = param(7); Cdc = param(8); Pmpp = param(9);

    DELTA = asin(2*Xg*Pmpp/(3*V0^2));
    PF = Pmpp;
    VDC = mp/kp*(Pmpp - Pref) + vdcref;
    Gpd = (3*V0^2*cos(DELTA))/(2*Xg);

    P = x(1);
    vdc = x(2);

    dx = zeros(2,1);
    dx(1) = -Gpd*mp*(P - PF) + Gpd*kp*(vdc - VDC);
    dx(2) = -(P - Pmpp)/(Cdc*VDC);
end
