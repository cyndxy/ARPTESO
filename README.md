clear; clc; close all;
fprintf('Starting simulation...\n');
t_span = [0 50];          
M = 4;                   
N_total = 8;              
I3 = eye(3);              

% --- 2. 通信拓扑
L1 = [ 2 -1  0  0;
      -1  3 -1  0;
       0 -1  3 -1;
       0  0 -1  2];
L2 = [-1  0  0  0;
       0 -1  0  0;
       0  0 -1  0;
       0  0  0 -1];
L1_inv_L2 = inv(L1) * L2; 
eta1_0 = [3.145; -0.55; 0.45];   
eta2_0 = [2.68; -0.2; 0.55];      
eta3_0 = [1.84; -0.2; pi/4];     
eta4_0 = [1.47; -0.49; 0.5];     
v1_0 = [0.10; 0.00; 0.00];
v2_0 = [0.08; 0.00; 0.00];
v3_0 = [0.06; 0.00; 0.00];
v4_0 = [0.08; 0.00; 0.00];
w1_0 = R_psi(eta1_0(3)) * v1_0;
w2_0 = R_psi(eta2_0(3)) * v2_0;
w3_0 = R_psi(eta3_0(3)) * v3_0;
w4_0 = R_psi(eta4_0(3)) * v4_0;
Lambda0 = 8*ones(3,1);   % RPTESO自适应律初值，满足 Lambda(0) >= 1

z0_1 = [eta1_0; w1_0; eta1_0; zeros(3,1); zeros(3,1); Lambda0];
z0_2 = [eta2_0; w2_0; eta2_0; zeros(3,1); zeros(3,1); Lambda0];
z0_3 = [eta3_0; w3_0; eta3_0; zeros(3,1); zeros(3,1); Lambda0];
z0_4 = [eta4_0; w4_0; eta4_0; zeros(3,1); zeros(3,1); Lambda0];

% 系统总初始状态
X0 = [z0_1; z0_2; z0_3; z0_4];

% --- 4. ODE求解器
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'MaxStep', 0.01);
[t, X] = ode45(@(t,X) ode_function(t, X, M, L1_inv_L2), t_span, X0, options);

fprintf('Simulation finished. Processing and plotting results...\n');
num_states_per_ship = 18;
eta_hist = zeros(3, M, length(t));   
e_hist = zeros(3, M, length(t));     
w_err_hist = zeros(3, M, length(t)); 
f_err_hist = zeros(3, M, length(t)); 
tau_hist = zeros(3, M, length(t));   
lambda_ppc_hist = zeros(3, length(t));
Lambda_obs_hist = zeros(3, M, length(t));

for k = 1:length(t)
    [~, plot_data] = ode_function(t(k), X(k,:)', M, L1_inv_L2);
    for i = 1:M
        eta_hist(:,i,k) = plot_data(i).eta_i;
        e_hist(:,i,k) = plot_data(i).e_i;
        w_err_hist(:,i,k) = plot_data(i).w_err;
        f_err_hist(:,i,k) = plot_data(i).f_err;
        tau_hist(:,i,k) = plot_data(i).tau_i;
        Lambda_obs_hist(:,i,k) = plot_data(i).Lambda_obs;
    end
    lambda_ppc_hist(:,k) = plot_data(1).lambda_t;
end

% ---------------------- IAE/ISE/ITAE指标计算 ----------------------
IAE = zeros(3, M); ISE = zeros(3, M); ITAE = zeros(3, M);
for j = 1:3
    for i = 1:M
        e_ij = squeeze(e_hist(j,i,:));
        IAE(j,i)  = trapz(t, abs(e_ij));
        ISE(j,i)  = trapz(t, e_ij.^2);
        ITAE(j,i) = trapz(t, t.*abs(e_ij));
    end
end

fprintf('\n==================== 误差性能指标 ====================\n');
err_names = {'x方向位置误差 (m)', 'y方向位置误差 (m)', '艏向角误差 (rad)'};
for j = 1:3
    fprintf('\n--- %s ---\n', err_names{j});
    fprintf('跟随船\t\tIAE\t\tISE\t\tITAE\n');
    for i = 1:M
        fprintf('Follower %d\t%.4f\t%.4f\t%.4f\n', i, IAE(j,i), ISE(j,i), ITAE(j,i));
    end
    fprintf('均值\t\t%.4f\t%.4f\t%.4f\n', mean(IAE(j,:)), mean(ISE(j,:)), mean(ITAE(j,:)));
end
fprintf('\n--- 全局均值（所有维度+所有跟随船） ---\n');
fprintf('IAE均值:  %.4f\n', mean(IAE(:)));
fprintf('ISE均值:  %.4f\n', mean(ISE(:)));
fprintf('ITAE均值: %.4f\n', mean(ITAE(:)));
fprintf('======================================================\n\n');


% ===================== 编队轨迹图 =====================
figure('Name', '编队轨迹', 'Position', [100,100,800,600]);
hold on; box on; axis equal;

leader_colors_traj = [
    141/255, 155/255, 206/255;   % 1
    112/255, 199/255, 212/255;   % 2
    123/255, 170/255, 203/255;   % 3
    184/255, 216/255, 150/255    % 4
];
colors = [
    247/255, 180/255, 101/255;   % 5
     193/255, 228/255, 232/255;     % 6
    247/255, 196/255, 191/255;   % 7
    234/255, 127/255, 109/255% 8
];

% 领航船标记形状
leader_marks = {'p', 's', 'd', '^'};

% 1. 绘制领航船轨迹（虚线）
t_plot = linspace(0, t(end), 200);
eta_l_path = get_leader_trajectories(t_plot);   
for l = 1:4
    plot(squeeze(eta_l_path(1,l,:)), squeeze(eta_l_path(2,l,:)), ...
         'Color', leader_colors_traj(l,:), 'LineStyle', '--', 'LineWidth', 1.2, ...
         'DisplayName', ['Leader' num2str(l)]);
end

% 2. 绘制跟随船轨迹（实线）
for i = 1:M
    plot(squeeze(eta_hist(1,i,:)), squeeze(eta_hist(2,i,:)), ...
         'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', ['USV' num2str(i)]);
end

time_points = [5, 20, 35, 50];  
leader_dot_colors = leader_colors_traj;
follower_dot_colors = colors;

for k = 1:length(time_points)
    tp = time_points(k);
    [~, idx] = min(abs(t - tp));
    eta_leader_now = get_leader_trajectories(t(idx));   
    leader_x = squeeze(eta_leader_now(1,:));
    leader_y = squeeze(eta_leader_now(2,:));
    K = convhull(leader_x, leader_y);
    plot(leader_x(K), leader_y(K), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    for l = 1:4
        plot(leader_x(l), leader_y(l), 'o', ...
             'MarkerSize', 12, 'MarkerFaceColor', leader_dot_colors(l,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    for i = 1:M
        plot(eta_hist(1,i,idx), eta_hist(2,i,idx), 'o', ...
             'MarkerSize', 12, 'MarkerFaceColor', follower_dot_colors(i,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    center_x = mean(leader_x);
    center_y = mean(leader_y);
    text(center_x, center_y, ['  t=', num2str(tp), 's'], ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', ...
         'EdgeColor', 'k', 'Margin', 2, 'HandleVisibility', 'off');
end

xlabel('x (m)', 'Interpreter', 'none'); ylabel('y (m)', 'Interpreter', 'none'); ylim([-4.5,4.5]);xlim([-4.5,4.5]);
legend_str = {'Leader1', 'Leader2', 'Leader3', 'Leader4', 'USV1', 'USV2', 'USV3', 'USV4'};
hLegend = legend(legend_str, 'Location', 'southoutside', 'Orientation', 'horizontal');
set(hLegend, 'Interpreter', 'latex');
hold off;
%包容跟踪误差
err_names = {'x方向位置误差', 'y方向位置误差', '艏向角误差'};
ylabels_latex = {'${e}_{x}$', '${e}_{y}$', '${e}_{\psi}$'};  
figure('Name', '包容跟踪误差', 'Position', [200, 200, 800, 900]);
line_handles = [];

for j = 1:3
    subplot(3,1,j);
    hold on; box on;
    h_up = plot(t, lambda_ppc_hist(j,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'upper bound');
    h_low = plot(t, -lambda_ppc_hist(j,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'lower bound');
    for i = 1:M
        plot(t, squeeze(e_hist(j,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', ['USV ' num2str(i)]);
    end
    if j == 1
        line_handles = [h_up, h_low];
        for i = 1:M
        end
    end
    
    xlabel('Time (s)', 'Interpreter', 'none');
    ylabel(ylabels_latex{j}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    xlim([0, 10]);
end


h_ax = subplot(3,1,1);  

all_lines = findobj(h_ax, 'Type', 'line');

usv_lines = [];
for i = 1:M
    usv_lines = [usv_lines, findobj(h_ax, 'DisplayName', ['USV ' num2str(i)])];
end
legend_handles = [findobj(h_ax, 'DisplayName', 'upper bound'), ...
                  findobj(h_ax, 'DisplayName', 'lower bound'), ...
                  usv_lines];

legend_str = {'upper bound', 'lower bound', 'USV1', 'USV2', 'USV3', 'USV4'};


hLegend = legend(legend_handles, legend_str, ...
                 'Location', 'southoutside', 'Orientation', 'horizontal');
set(hLegend, 'Interpreter', 'latex', 'FontSize', 10);

%增加原始指数型性能函数对比）
err_names = {'x方向位置误差', 'y方向位置误差', '艏向角误差'};
ylabels_latex = {'${e}_{x}$', '${e}_{y}$', '${e}_{\psi}$'};
figure('Name', '包容跟踪误差', 'Position', [200, 200, 800, 900]);


lambda0_exp = [1.0; 1.0; 1.0];
lambda_inf_exp = [0.1; 0.1; 0.1];
beta_exp = 1;

exp_lambda = (lambda0_exp - lambda_inf_exp) * exp(-beta_exp * t)' + lambda_inf_exp;  % 3×N
exp_lambda_neg = -exp_lambda;

for j = 1:3
    subplot(3,1,j);
    hold on; box on;
  
    h_up_imp = plot(t, lambda_ppc_hist(j,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'upper bound');
    h_low_imp = plot(t, -lambda_ppc_hist(j,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'lower bound');
   
    h_up_exp = plot(t, exp_lambda(j,:), 'k-', 'LineWidth', 1.2, 'DisplayName', 'upper bound');
    h_low_exp = plot(t, exp_lambda_neg(j,:), 'k--', 'LineWidth', 1.2, 'DisplayName', 'lower');
   
    for i = 1:M
        plot(t, squeeze(e_hist(j,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', ['USV ' num2str(i)]);
    end
    
    xlabel('Time (s)', 'Interpreter', 'none');
    ylabel(ylabels_latex{j}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'latex');
    xlim([0, 10]);
end


h_ax = subplot(3,1,1);

usv_lines = [];
for i = 1:M
    usv_lines = [usv_lines, findobj(h_ax, 'DisplayName', ['USV ' num2str(i)])];
end
h_up_imp = findobj(h_ax, 'DisplayName', 'upper bound');
h_low_imp = findobj(h_ax, 'DisplayName', 'upper bound');
h_up_exp = findobj(h_ax, 'DisplayName', 'lower bound');
h_low_exp = findobj(h_ax, 'DisplayName', 'lower bound');

legend_handles = [h_up_imp(:)', h_low_imp(:)', h_up_exp(:)', h_low_exp(:)', usv_lines(:)'];
legend_str = {'upper bound', 'lower bound', 'upper bound', 'lower bound', ...
              'USV1', 'USV2', 'USV3', 'USV4'};

hLegend = legend(legend_handles, legend_str, ...
                 'Location', 'southoutside', 'Orientation', 'horizontal');
set(hLegend, 'Interpreter', 'latex', 'FontSize', 10);

% ===================== 合并后的 3x2 观测误差图（速度 + 干扰） =====================
figure('Name', '观测误差（速度与干扰）', 'Position', [300, 200, 1000, 900]);

% 子图标签
vel_labels = {'$\tilde{z}_{2u}$', '$\tilde{z}_{2v}$', '$\tilde{z}_{2\psi}$'};
dist_labels = {'$\tilde{z}_{3u}$', '$\tilde{z}_{3v}$', '$\tilde{z}_{3\psi}$'};

for j = 1:3
    % ---------- 左列：速度观测误差 ----------
    subplot(3,2,2*j-1);
    hold on; box on;
    for i = 1:M
        plot(t, squeeze(w_err_hist(j,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', ['USV' num2str(i)]);
    end
    xline(0.8, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{T}_{\mathcal{P}}=0.8$');
    xlabel('Time (s)', 'Interpreter', 'none');
    ylabel(vel_labels{j}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'latex');
    xlim([0, 10]);

    % ---------- 右列：干扰观测误差 ----------
    subplot(3,2,2*j);
    hold on; box on;
    for i = 1:M
        plot(t, squeeze(f_err_hist(j,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', ['USV' num2str(i)]);
    end
    xline(0.8, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{T}_{\mathcal{P}}=0.8$');
    xlabel('Time (s)', 'Interpreter', 'none');
    ylabel(dist_labels{j}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'latex');
    xlim([0, 10]);
end

ax_first = subplot(3,2,1);
usv_lines = [];
for i = 1:M
    usv_lines = [usv_lines, findobj(ax_first, 'DisplayName', ['USV' num2str(i)])];
end
tp_line = findobj(ax_first, 'DisplayName', '$\mathcal{T}_{\mathcal{P}}=0.8$');
legend_handles = [usv_lines, tp_line];
legend_str = {'USV1', 'USV2', 'USV3', 'USV4', '$\mathcal{T}_{\mathcal{P}}=0.8$'};
hLegend = legend(legend_handles, legend_str, ...
                 'Location', 'southoutside', 'Orientation', 'horizontal');
set(hLegend, 'Interpreter', 'latex', 'FontSize', 10);

figure('Name', 'RPTESO自适应增益', 'Position', [600, 250, 1200, 400]);

% 预定义各子图的y轴标签
gamma_labels = {'Adaptive gains ($x$)', 'Adaptive gains ($y$)', 'Adaptive gains ($\psi$)'};

for j = 1:3
    subplot(1,3,j);
    hold on; box on;
    for i = 1:M
        plot(t, squeeze(Lambda_obs_hist(j,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', ['USV' num2str(i)]);
    end
    xlabel('Time (s)', 'Interpreter', 'none');
    ylabel(gamma_labels{j}, 'Interpreter', 'latex', ...
           'FontSize', 10, 'FontWeight', 'bold');
    xlim([0, 4]);
end
subplot(1,3,1);
usv_lines_gamma = [];
for i = 1:M
    usv_lines_gamma = [usv_lines_gamma, findobj(gca, 'DisplayName', ['USV' num2str(i)])];
end
legend_str_gamma = {'USV1', 'USV2', 'USV3', 'USV4'};
hLegend_gamma = legend(usv_lines_gamma, legend_str_gamma, ...
                       'Location', 'southoutside', 'Orientation', 'horizontal');
set(hLegend_gamma, 'Interpreter', 'latex', 'FontSize', 10);

fprintf('Plotting complete.\n');

% ===================== 函数定义 =====================

function [M, C_v, D_v] = get_ship_dynamics(v)
    % Cybership II 动力学模型（保持不变）
    m = 23.8; xg = 0.046; Iz = 1.760;
    X_udot = -2.0; Y_vdot = -10.0; Y_rdot = 0.0; N_vdot = 0.0; N_rdot = -1.0;
    M_RB = [m, 0, m*xg; 0, m, 0; m*xg, 0, Iz];
    M_A = [X_udot, 0, 0; 0, Y_vdot, Y_rdot; 0, N_vdot, N_rdot];
    M = M_RB + M_A;

    u = v(1); v_sway = v(2); r = v(3);
    C_RB = [0, 0, -m*(xg*r+v_sway); 0, 0, m*u; m*(xg*r+v_sway), -m*u, 0];
    C_A = [0, 0, -Y_vdot*v_sway - (Y_rdot+N_vdot)*r; 0, 0, X_udot*u; ...
           Y_vdot*v_sway+(Y_rdot+N_vdot)*r, -X_udot*u, 0];
    C = C_RB + C_A;
    C_v = C * v;

    Xu = -0.72253; Yv = -0.88965; Yr = -7.250; Nv = 0.03130; Nr = -1.900;
    D_lin = [Xu, 0, 0; 0, Yv, Yr; 0, Nv, Nr];
    X_uu = -1.32742; Y_vv = -36.47287; Y_rr = -3.450; Y_rv = -0.805; Y_vr = -0.845;
    N_vv = 3.95645; N_rr = -0.750; N_rv = 0.130;
    D_nlin = [X_uu*abs(u), 0, 0; ...
              0, Y_vv*abs(v_sway)+Y_rv*abs(r), Y_rr*abs(r)+Y_vr*abs(v_sway); ...
              0, N_vv*abs(v_sway)+N_rv*abs(r), N_rr*abs(r)];
    D = D_lin + D_nlin;
    D_v = D * v;
end

function y = sig(x, alpha)
    delta = 1e-12;
    y = sign(x) .* (abs(x) + delta).^alpha;
end

function y = dsig(x, alpha)
    delta = 1e-12;
    y = alpha .* (abs(x) + delta).^(alpha - 1);
end

function y = sat_vec(x)
    y = max(-1, min(1, x));
end

function [eta_leaders, w_leaders, w_dot_leaders] = get_leader_trajectories(t)
    w1 = pi/30;
    N = length(t);
    num_leaders = 4;
    eta_leaders = zeros(3, num_leaders, N);
    w_leaders = zeros(3, num_leaders, N);
    w_dot_leaders = zeros(3, num_leaders, N);

    for k = 1:N
        tk = t(k);
        % 领航船1
        r1 = 3.56; phi1 = -0.322;
        eta_leaders(:,1,k) = [r1*cos(w1*tk+phi1); r1*sin(w1*tk+phi1); 0];
        w_leaders(:,1,k) = [-r1*w1*sin(w1*tk+phi1); r1*w1*cos(w1*tk+phi1); 0];
        w_dot_leaders(:,1,k) = [-r1*w1^2*cos(w1*tk+phi1); -r1*w1^2*sin(w1*tk+phi1); 0];

        % 领航船2
        r2 = 3.8; phi2 = 0;
        eta_leaders(:,2,k) = [r2*cos(w1*tk+phi2); r2*sin(w1*tk+phi2); 0];
        w_leaders(:,2,k) = [-r2*w1*sin(w1*tk+phi2); r2*w1*cos(w1*tk+phi2); 0];
        w_dot_leaders(:,2,k) = [-r2*w1^2*cos(w1*tk+phi2); -r2*w1^2*sin(w1*tk+phi2); 0];

        % 领航船3
        r3 = 1.5; phi3 = 0;
        eta_leaders(:,3,k) = [r3*cos(w1*tk+phi3); r3*sin(w1*tk+phi3); 0];
        w_leaders(:,3,k) = [-r3*w1*sin(w1*tk+phi3); r3*w1*cos(w1*tk+phi3); 0];
        w_dot_leaders(:,3,k) = [-r3*w1^2*cos(w1*tk+phi3); -r3*w1^2*sin(w1*tk+phi3); 0];

        % 领航船4
        r4 = 1.31; phi4 = -0.7096;
        eta_leaders(:,4,k) = [r4*cos(w1*tk+phi4); r4*sin(w1*tk+phi4); 0];
        w_leaders(:,4,k) = [-r4*w1*sin(w1*tk+phi4); r4*w1*cos(w1*tk+phi4); 0];
        w_dot_leaders(:,4,k) = [-r4*w1^2*cos(w1*tk+phi4); -r4*w1^2*sin(w1*tk+phi4); 0];
    end
end

function [dXdt, plot_data] = ode_function(t, X, M_followers, L1_inv_L2)
    I3 = eye(3);
    num_states_per_ship = 18;
    dXdt = zeros(size(X));
    plot_data = struct('eta_i', [], 'e_i', [], 'w_err', [], 'f_err', [], ...
                       'tau_i', [], 'lambda_t', [], 'Lambda_obs', []);

    [eta_leaders, w_leaders, w_dot_leaders] = get_leader_trajectories(t);
    eta_id_vec = reshape(eta_leaders, [], 1);
    w_leaders_vec = reshape(w_leaders, [], 1);
    w_dot_id_vec = reshape(w_dot_leaders, [], 1);

    % ====================== PTSM 参数 ======================
    alpha_s = 0.5;
    beta_s = 1;
    T_s = 1;
    alpha_u = 0.5;
    beta_u = 1;
    T_u = 1;
    epsilon0 = 0.1;
    l1 = 2 * epsilon0^(-alpha_s);
    l2 = -epsilon0^(-2*alpha_s);
    bar_alpha_s = pi / (2^(1-alpha_s/2) * alpha_s * beta_s * T_s);
    bar_beta_s = beta_s * pi / (2^(1+alpha_s/2) * alpha_s * T_s);
    bar_alpha_u = pi / (2^(1-alpha_u/2) * alpha_u * beta_u * T_u);
    bar_beta_u = beta_u * pi / (2^(1+alpha_u/2) * alpha_u * T_u);
    % ======================================================

    for i = 1:M_followers
        idx_start = (i-1)*num_states_per_ship + 1;
        idx_end = i*num_states_per_ship;
        X_i = X(idx_start:idx_end);

        eta_i = X_i(1:3);
        w_i = X_i(4:6);
        eta_hat_i = X_i(7:9);
        w_hat_i = X_i(10:12);
        f_hat_i = X_i(13:15);
        Lambda_obs_i = X_i(16:18);

        psi_i = eta_i(3);
        R_i = R_psi(psi_i);
        v_i = R_i' * w_i;

        [M_ship, C_ship_v, D_ship_v] = get_ship_dynamics(v_i);

        % 环境干扰（保持不变）
        d_i = [1.3 + 2*sin(0.02*t) + 1.5*sin(0.1*t);
               0.9 + 2*sin(0.02*t - pi/6) + 1.5*sin(0.3*t);
               -2*sin(0.09*t + pi/3) - 4*sin(0.01*t)];

        R_dot_i = [-sin(psi_i)*w_i(3), -cos(psi_i)*w_i(3), 0;
                    cos(psi_i)*w_i(3), -sin(psi_i)*w_i(3), 0;
                    0, 0, 0];
        f_i = R_dot_i * v_i + R_i * (M_ship \ (-C_ship_v - D_ship_v + d_i));

        % ====================== 【仅修改此处：改进性能函数】 ======================
        lambda_0   = [1.0; 1.0; 1.0];
        lambda_inf = [0.1; 0.1; 0.1];

        T_p = 2.0;

        if t < T_p
            c = cos(pi*t/T_p);
            s = sin(pi*t/T_p);
            phi = (1 + c)/2;

            lambda_t      = (lambda_0 - lambda_inf) .* (phi.^2) + lambda_inf;
            lambda_dot_t  = -(lambda_0 - lambda_inf) .* (pi/T_p) .* phi .* s;
            lambda_ddot_t = -(lambda_0 - lambda_inf) .* (pi^2/(2*T_p^2)) .* (c + cos(2*pi*t/T_p));
        else
            lambda_t      = lambda_inf;
            lambda_dot_t  = zeros(3,1);
            lambda_ddot_t = zeros(3,1);
        end
        % =======================================================================

        e_i = eta_i + kron(L1_inv_L2(i,:), I3) * eta_id_vec;
        e_dot_i = w_i + kron(L1_inv_L2(i,:), I3) * w_leaders_vec;

        delta_ppc = 1e-9;
        e_i_safe = max(-lambda_t + delta_ppc, min(lambda_t - delta_ppc, e_i));

        eps_i = 0.5 * (log(lambda_t + e_i_safe) - log(lambda_t - e_i_safe));
        eps_dot_i = (lambda_t.*e_dot_i - e_i_safe.*lambda_dot_t) ./ (lambda_t.^2 - e_i_safe.^2);

        % 计算 B1, B2, B3 
        tanh_eps = e_i_safe ./ lambda_t;
        sech2_eps = (lambda_t.^2 - e_i_safe.^2) ./ lambda_t.^2;
        B1_vec = lambda_ddot_t .* tanh_eps + 2*lambda_dot_t .* sech2_eps .* eps_dot_i;
        B2_vec = 2*lambda_t .* tanh_eps .* sech2_eps .* (eps_dot_i.^2);
        B3_inv = diag(lambda_t ./ (lambda_t.^2 - e_i_safe.^2));
        B3 = diag((lambda_t.^2 - e_i_safe.^2) ./ lambda_t);

        % ====================== PTSM 滑模面计算 ======================
        abs_eps = abs(eps_i);
        f_eps = zeros(3,1);
        dot_f_eps = zeros(3,1);
        for j = 1:3
            if abs_eps(j) >= epsilon0
                f_eps(j) = sig(eps_i(j), 1-alpha_s);
                dot_f_eps(j) = (1-alpha_s) * sig(eps_i(j), -alpha_s) * eps_dot_i(j);
            else
                f_eps(j) = l1 * eps_i(j) + l2 * sig(eps_i(j), 1+alpha_s);
                dot_f_eps(j) = l1 * eps_dot_i(j) + l2 * (1+alpha_s) * sig(eps_i(j), alpha_s) * eps_dot_i(j);
            end
        end

        s_i = eps_dot_i + bar_alpha_s * f_eps + bar_beta_s * sig(eps_i, 1+alpha_s);
        % ============================================================

        % ====================== 预定义时间趋近律 ======================
        delta_s = 0.01;
        abs_s = abs(s_i);
        sig_s1 = (abs_s < delta_s) .* (s_i ./ (delta_s.^(1-alpha_u))) + ...
                 (abs_s >= delta_s) .* sig(s_i, 1-alpha_u);
        sig_s2 = (abs_s < delta_s) .* (s_i ./ (delta_s.^(1+alpha_u))) + ...
                 (abs_s >= delta_s) .* sig(s_i, 1+alpha_u);
        term_slide = bar_alpha_u * sig_s1 + bar_beta_u * sig_s2;
        % ============================================================

        % ====================== 控制律 ======================
        term_eps = bar_alpha_s * dot_f_eps + bar_beta_s * (1+alpha_s) * (sig(eps_i, alpha_s) .* eps_dot_i);
        e_ddot_ref = kron(L1_inv_L2(i,:), I3) * w_dot_id_vec;
        tau_i = M_ship * R_i' * ( -B3 * ( term_slide + (B1_vec + B2_vec) + term_eps ) + e_ddot_ref - f_hat_i );
        % ====================================================

        % ====================== RPTESO 观测器======================
        z_tilde_i1 = eta_i - eta_hat_i;

        g1_obs = 6; g2_obs = 8; g3_obs = 2;
        K_obs  = 4; k3_obs = 4;
        k1_obs = 3.34 * k3_obs^(1/3);
        k2_obs = 5.3  * k3_obs^(2/3);
        Theta_obs = 1e2 * ones(3,1);
        tf_obs = 0.8;
        Xi_max = 100;

        if t < tf_obs
            l_sw = 1;
            Xi_t = sec(pi*t/(2*tf_obs))^2;
            Xi_t = min(Xi_t, Xi_max);
        else
            l_sw = 0;
            Xi_t = 1;
        end

        d_Lambda_obs_i = abs(z_tilde_i1 - Theta_obs .* sat_vec(Theta_obs .* z_tilde_i1));
        Lambda_safe = max(Lambda_obs_i, 1e-6);

        term1_rpteso = l_sw * ( Lambda_safe .* Xi_t .* g1_obs .* z_tilde_i1 ...
                      + (d_Lambda_obs_i ./ Lambda_safe) .* z_tilde_i1 );
        term2_rpteso = l_sw * ( (Lambda_safe.^2) .* (Xi_t^2) .* g2_obs .* z_tilde_i1 );
        term3_rpteso = l_sw * ( (Lambda_safe.^3) .* (Xi_t^3) .* g3_obs .* z_tilde_i1 );

        term1_hosm = (1 - l_sw) * k1_obs * K_obs^(1/3) * sig(z_tilde_i1, 2/3);
        term2_hosm = (1 - l_sw) * k2_obs * K_obs^(2/3) * sig(z_tilde_i1, 1/3);
        term3_hosm = (1 - l_sw) * k3_obs * K_obs * sign(z_tilde_i1);

        d_eta_hat_i = w_hat_i + term1_rpteso + term1_hosm;
        d_w_hat_i = R_i*(M_ship \ tau_i) + f_hat_i + term2_rpteso + term2_hosm;
        d_f_hat_i = term3_rpteso + term3_hosm;
   
% ================================================================
        d_eta_i = w_i;
        d_w_i = R_i*(M_ship \ tau_i) + f_i;

        dXdt(idx_start:idx_end) = [d_eta_i; d_w_i; d_eta_hat_i; d_w_hat_i; d_f_hat_i; d_Lambda_obs_i];

        plot_data(i).eta_i = eta_i;
        plot_data(i).e_i = e_i;
        plot_data(i).w_err = w_i - w_hat_i;
        plot_data(i).f_err = f_i - f_hat_i;
        plot_data(i).tau_i = tau_i;
        plot_data(i).lambda_t = lambda_t;
        plot_data(i).Lambda_obs = Lambda_obs_i;
    end
end

function R = R_psi(psi)
    R = [cos(psi) -sin(psi) 0;
         sin(psi)  cos(psi) 0;
         0         0        1];
end
