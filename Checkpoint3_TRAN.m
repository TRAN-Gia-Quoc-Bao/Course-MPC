%% MPC checkpoint 3 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
% Settings
Nsim = 30;                                                      % number of samples we simulate
Np = 2;                                                         % prediction horizon
lesx = zeros(Nsim, 1); lesx(1) = 0;                             % state initialization
lesu = zeros(Nsim, 1);                                          % control initialization
umin = -0.4; umax = +0.6;                                       % control bounds
xd = 2;                                                         % reference
alphaValues = [0 1 10];                                         % loop with different values for alpha
legendsStrings = cell(length(alphaValues), 1);
for j = 1 : length(alphaValues)                                                  
    alpha = alphaValues(j);                                                 

    % MPC
    for i = 1 : Nsim - 1
        lesu(i) = Kdex(lesx(i), Np, umin, umax, xd, alpha);     % calculate control using optimization
        lesx(i+1) = lesx(i) + lesu(i);                          % x(k + 1) = x(k) + u(k)                  
    end
    lesu(Nsim) = Kdex(lesx(Nsim), Np, umin, umax, xd, alpha);   % the last sample

    subplot(211);
    plot(lesx, '-o', 'LineWidth', 2);
    grid on;
    hold on;
    xlim([1 Nsim]);
    xlabel('Sample', 'FontSize', 16);
    ylabel('State', 'FontSize', 16);
    title('State evolution', 'FontSize', 16);
    subplot(212);
    stairs(lesu, '-o', 'LineWidth', 2);
    grid on;
    hold on;
    xlim([1 Nsim]);
    xlabel('Sample', 'FontSize', 16);
    ylabel('Control', 'FontSize', 16);
    title('Control evolution', 'FontSize', 16);
    legendsStrings{j} = ['alpha = ', num2str(alpha)];
end
legend(legendsStrings, 'Interpreter', 'none');

%% Functions
function u = Kdex(x, Np, umin, umax, xd, alpha)
    inter = p_unc(x, Np, xd, alpha);                            % optimization
    u = max(umin, min(umax, inter));                            % apply the bounds
    return
end
%--------------------------------------------
function f = p_unc(x, Np, xd, alpha)
    inter1 = psi1(Np);
    inter0 = psi0(x, Np, xd);
    f = -inter1'*inter0/(inter1'*inter1 + alpha);               % answer of the optimization problem
    return
end
%--------------------------------------------
function f = psi1(Np)
    f = ones(Np, 1);
    return
end
%-------------------------------------------
function f = psi0(x, Np, xd)
    f = ones(Np, 1)*(x - xd);
    return
end
%-------------------------------------------