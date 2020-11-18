%% MPC checkpoint 8 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
% System
tau = 0.25; N = 20; xd = [2; 0; 0]; ud = 0; Q = eye(3); R = 1; P = 10*eye(3);
sysC = ss([0 1 0; 0 0 1; 0 0 0], [0; 0; 1], [1 0 0], [0]);
sysD = c2d(sysC, tau);
% Constraints
state2ConstraintValues = [0.3 0.25 0.15];                                   % loop with different values for constraint
legendsStrings = cell(length(state2ConstraintValues), 1);
for j = 1 : length(state2ConstraintValues)                                                                                                   
    umax = [0.45]; umin = -umax;
    xmax = [inf; state2ConstraintValues(j); 0.25]; xmin = -xmax; 
    [Aineq, Bineq] = construct_constraints_mat(umax, umin, xmax, xmin, N);  % if varying constraints: put this inside loop
    [Phi, Psi] = construct_Gamma(sysD.A, sysD.B, N);

    % Cost function
    [H, h, c] = compute_cost_matrices(Q, R, P, N, xd, ud);

    % The open-looped case (cannot go beyond the prediction horizon)
%     x0 = [0; 0; 0];
%     z = quadprog(2*H, h, Aineq, Bineq, Psi, -Phi*x0);
%     u = []; out = [];
%     for k = 0 : N - 1
%         u = [u; z(4*k + 1)];
%         out = [out; z(4*k + 2)];
%     end
%     figure('Name', 'Open-loop system');
%     subplot(211);
%     plot(out, '-o', 'LineWidth', 2);
%     grid on;
%     hold on;
%     xlim([1 N]);
%     xlabel('Sample', 'FontSize', 16);
%     ylabel('Output', 'FontSize', 16);
%     title('Open-loop system: output evolution', 'FontSize', 16);
%     subplot(212);
%     stairs(u, '-o', 'LineWidth', 2);
%     grid on;
%     hold on;
%     xlim([1 N]);
%     xlabel('Sample', 'FontSize', 16);
%     ylabel('Control', 'FontSize', 16);
%     title('Open-loop system: control evolution', 'FontSize', 16);

    % The closed-looped case (can go to infinity)
    Nsim = 100;                                                             % number of samples we simulate
    lesx = zeros(Nsim, 3); lesx(1, :) = 0;                                  % state initialization
    lesu = zeros(Nsim, 1);                                                  % control initialization
    for i = 1 : Nsim - 1
        z = quadprog(2*H, h, Aineq, Bineq, Psi, -Phi*lesx(i, :)');          % the current state is in the equality constraint
        lesu(i) = z(1);                                                     % use only the first part of the control
        lesx(i + 1, :) = (sysD.A*lesx(i, :)' + sysD.B*lesu(i))';            % x(k + 1) = A*x(k) + B*u(k)                  
    end
    %figure('Name', 'Closed-loop system');
    subplot(211);
    plot(lesx(:, 1), '-o', 'LineWidth', 2);
    grid on;
    hold on;
    xlim([1 Nsim]);
    xlabel('Sample', 'FontSize', 16);
    ylabel('Output', 'FontSize', 16);
    title('Closed-loop system: output evolution', 'FontSize', 16);
    subplot(212);
    stairs(lesu, '-o', 'LineWidth', 2);
    grid on;
    hold on;
    xlim([1 Nsim]);
    xlabel('Sample', 'FontSize', 16);
    ylabel('Closed-loop system: control', 'FontSize', 16);
    title('Closed-loop system: control evolution', 'FontSize', 16);
    legendsStrings{j} = ['Constraint on state 2 = ', num2str(state2ConstraintValues(j))];
end
legend(legendsStrings, 'Interpreter', 'none');

%% Functions
function [Aineq, Bineq] = construct_constraints_mat(umax, umin, xmax, xmin, N)
    nx = length(xmax);
    nu = length(umax);
    Aineq = [];
    Bineq = [];
    for k = 0 : N - 1
        sub = zeros(1, N);                                                  % initiate the sub matrix
        sub(k + 1) = 1;                                                     % move 1 along the sub matrix
        Sk = kron(sub, eye(nx + nu));                                       % calculate Sk
        Gk = [eye(nu) zeros(nu, nx); -eye(nu) zeros(nu, nx); zeros(nx, nu) eye(nx); zeros(nx, nu) -eye(nx)];
        dk = [umax; -umin; xmax; -xmin];
        Aineq = [Aineq; Gk*Sk];                                             % stock the matrix Aineq
        Bineq = [Bineq; dk];                                                % stock the matrix Bineq
    end
    return
end
%--------------------------------------------
function [Phi, Psi] = construct_Gamma(A, B, N)
    nx = size(A, 1);
    nu = size(B, 2);
    lineJ = zeros(nx, nx + N*(nx + nu));
    Gam = [];
    for j = 0 : N - 1
        lineJ(1 : nx, 1 + j*(nx + nu) : nx + j*(nx + nu)) = A; 
        lineJ(1 : nx, 1 + nx + j*(nx + nu) : nx + nu + j*(nx + nu)) = B; 
        lineJ(1 : nx, 1 + (nx + nu) + j*(nx + nu) : 2*nx + nu + j*(nx + nu)) = -eye(nx); 
        Gam = [Gam; lineJ];
        lineJ = zeros(nx, nx + N*(nx + nu));
    end
    Phi = Gam(:, 1 : nx);
    Psi = Gam(:, nx + 1 : end);
    return
end
%--------------------------------------------
function [H, h, c] = compute_cost_matrices(Q, R, P, N, xd, ud)
    nx = size(Q, 1);                                                        % for generality
    nu = size(R, 1);                                                        % for generality
    sub = zeros(1, N);                                                      % initiate the sub matrix 
    H = zeros(N*(nx + nu), N*(nx + nu));
    h = zeros(N*(nx + nu), 1);
    c = 0;
    for k = 0 : N - 2                                                       % loop until k = N - 2
        Qk = [R zeros(nu, nx); zeros(nx, nu) Q];
        hk = -2*[R*ud; Q*xd];
        ck = xd'*Q*xd + ud'*R*ud;
        sub(k + 1) = 1;                                                     % move 1 along the sub matrix
        Sk = kron(sub, eye(nx + nu));                                       % calculate Sk
        sub = zeros(1, N);                                                  % reset the sub matrix
        H = H + Sk'*Qk*Sk;
        h = h + (hk'*Sk)';
        c = c + ck;
    end
    % Add terminal penalty
    sub(N) = 1;                                                             % the last term is different
    Sk = kron(sub, eye(nx + nu));                                           % calculate the last Sk
    QN = [R zeros(nu, nx); zeros(nx, nu) Q + P];                            % the QN with added P for terminal penalty 
    hN = -2*[R*ud; (Q + P)*xd];                                             % the HN with added P for terminal penalty
    cN = xd'*(Q + P)*xd + ud'*R*ud;                                         % the CN with added P for terminal penalty
    H = H + Sk'*QN*Sk;
    h = h + (hN'*Sk)';
    c = c + cN;
    return
end