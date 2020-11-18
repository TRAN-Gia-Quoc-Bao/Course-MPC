%% MPC checkpoint 5 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
nx = 3;
nu = 2;
xd = [1; 2; 0];
ud = [-0.1; 0.1];
N = 2;
P = 10*eye(nx);
Q = [1; 2; 3]*[1 2 2] + diag([1; 0.2; 1.5]);
R = 0.15*eye(nu);
[H, h, c] = compute_cost_matrices(Q, R, P, N, xd, ud)

%% Functions
function [H, h, c] = compute_cost_matrices(Q, R, P, N, xd, ud)
    nx = size(Q, 1);                                    % for generality
    nu = size(R, 1);                                    % for generality
    sub = zeros(1, N);                                  % initiate the sub matrix 
    H = zeros(N*(nx + nu), N*(nx + nu));
    h = zeros(N*(nx + nu), 1);
    c = 0;
    for k = 0 : N - 2                                   % loop until k = N - 2
        Qk = [R zeros(nu, nx); zeros(nx, nu) Q];
        hk = -2*[R*ud; Q*xd];
        ck = xd'*Q*xd + ud'*R*ud;
        sub(k + 1) = 1;                                 % move 1 along the sub matrix
        Sk = kron(sub, eye(nx + nu));                   % calculate Sk
        sub = zeros(1, N);                              % reset the sub matrix
        H = H + Sk'*Qk*Sk;
        h = h + (hk'*Sk)';
        c = c + ck;
    end
    % Add terminal penalty
    sub(N) = 1;                                         % the last term is different
    Sk = kron(sub, eye(nx + nu));                       % calculate the last Sk
    QN = [R zeros(nu, nx); zeros(nx, nu) Q + P];        % the QN with added P for terminal penalty 
    hN = -2*[R*ud; (Q + P)*xd];                         % the HN with added P for terminal penalty
    cN = xd'*(Q + P)*xd + ud'*R*ud;                     % the CN with added P for terminal penalty
    H = H + Sk'*QN*Sk;
    h = h + (hN'*Sk)';
    c = c + cN;
    return
end