%% MPC checkpoint 7 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
umax = [3; 4]; umin = -[3; 4];
xmax = [5; 6; 7; 8]; xmin = -[5; 6; 7; 8]; 
N = 12;
[Aineq, Bineq] = construct_constraints_mat(umax, umin, xmax, xmin, N);

%% Functions
function [Aineq, Bineq] = construct_constraints_mat(umax, umin, xmax, xmin, N)
    nx = length(xmax);
    nu = length(umax);
    Aineq = [];
    Bineq = [];
    for k = 0 : N - 1
        sub = zeros(1, N);                                  % initiate the sub matrix
        sub(k + 1) = 1;                                     % move 1 along the sub matrix
        Sk = kron(sub, eye(nx + nu));                       % calculate Sk
        Gk = [eye(nu) zeros(nu, nx); -eye(nu) zeros(nu, nx); zeros(nx, nu) eye(nx); zeros(nx, nu) -eye(nx)];
        dk = [umax; -umin; xmax; -xmin];
        Aineq = [Aineq; Gk*Sk];                             % stock the matrix Aineq
        Bineq = [Bineq; dk];                                % stock the matrix Bineq
    end
    return
end