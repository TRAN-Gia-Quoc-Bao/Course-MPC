%% MPC checkpoint 6 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
tau = 0.25; N = 20; xd = [2; 0; 0]; ud = 0; Q = eye(3); R = 1; P = 10*eye(3);
% umax = [0.1]; umin = -umax;                                       % just to see this go wrong when we apply saturated unconstrained MPC instead of constrained MPC
sysC = ss([0 1 0; 0 0 1; 0 0 0], [0; 0; 1], [1 0 0], [0]);
sysD = c2d(sysC, tau);
Nsim = 50;                                                          % number of samples we simulate                                                        
lesx = zeros(Nsim, 3); lesx(1, :) = 0;                              % state initialization
lesu = zeros(Nsim, 1);                                              % control initialization
[K, v] = KMPC(sysD.A, sysD.B, N, Q, R, P, xd, ud);
for i = 1 : Nsim - 1
    lesu(i) = -K*lesx(i, :)' + v;                                   % use only the first part of the control
%     lesu(i) = max(umin, min(umax, lesu(i)));                        % just to see this go wrong when we apply saturated unconstrained MPC instead of constrained MPC
    lesx(i + 1, :) = (sysD.A*lesx(i, :)' + sysD.B*lesu(i))';        % x(k + 1) = A*x(k) + B*u(k)                  
end
lesu(Nsim, :) = -K*lesx(Nsim, :)' + v;                              % the last sample
subplot(211);
plot(lesx(:, 1), '-o', 'LineWidth', 2);
grid on;
hold on;
xlim([1 Nsim]);
xlabel('Sample', 'FontSize', 16);
ylabel('Output', 'FontSize', 16);
title('Unconstrained MPC: output evolution', 'FontSize', 16);
subplot(212);
stairs(lesu, '-o', 'LineWidth', 2);
grid on;
hold on;
xlim([1 Nsim]);
xlabel('Sample', 'FontSize', 16);
ylabel('Control', 'FontSize', 16);
title('Unconstrained MPC: control evolution', 'FontSize', 16);

%% Functions
% K_MPC
function [K, v] = KMPC(A, B, N, Q, R, P, xd, ud)
    [Phi, Psi] = construct_Gamma(A, B, N);
    [H, h, c] = compute_cost_matrices(Q, R, P, N, xd, ud);
    M = inv([2*H Psi'; Psi zeros(size(Psi, 1), size(Psi, 1))]);
    K = M(1, size(Psi, 2) + 1 : end)*Phi;
    v = -M(1, 1 : size(Psi, 2))*h;
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
    nx = size(Q, 1);                                                % for generality
    nu = size(R, 1);                                                % for generality
    sub = zeros(1, N);                                              % initiate the sub matrix 
    H = zeros(N*(nx + nu), N*(nx + nu));
    h = zeros(N*(nx + nu), 1);
    c = 0;
    for k = 0 : N - 2                                               % loop until k = N - 2
        Qk = [R zeros(nu, nx); zeros(nx, nu) Q];
        hk = -2*[R*ud; Q*xd];
        ck = xd'*Q*xd + ud'*R*ud;
        sub(k + 1) = 1;                                             % move 1 along the sub matrix
        Sk = kron(sub, eye(nx + nu));                               % calculate Sk
        sub = zeros(1, N);                                          % reset the sub matrix
        H = H + Sk'*Qk*Sk;
        h = h + (hk'*Sk)';
        c = c + ck;
    end
    % Add terminal penalty
    sub(N) = 1;                                                     % the last term is different
    Sk = kron(sub, eye(nx + nu));                                   % calculate the last Sk
    QN = [R zeros(nu, nx); zeros(nx, nu) Q + P];                    % the QN with added P for terminal penalty 
    hN = -2*[R*ud; (Q + P)*xd];                                     % the HN with added P for terminal penalty
    cN = xd'*(Q + P)*xd + ud'*R*ud;                                 % the CN with added P for terminal penalty
    H = H + Sk'*QN*Sk;
    h = h + (hN'*Sk)';
    c = c + cN;
    return
end