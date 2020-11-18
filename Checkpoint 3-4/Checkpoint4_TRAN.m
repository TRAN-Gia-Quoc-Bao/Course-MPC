%% MPC checkpoint 4 _ Gia Quoc Bao TRAN _ 2020

%% Default commands
close all; clear all; clc;

%% Main program
nx = 4;
nu = 2;
N = 3;
x0 = randn(nx, 1);                          % initial condition 
A = 0.3*randn(nx, nx);                      % A matrix
B = randn(nx, nu);                          % B matrix
[Phi, Psi] = construct_Gamma(A, B, N);      % Phi and Psi matrices
lesu = randn(N, nu);
xact = x0;
z = [];
for i = 1 : N
    uact = lesu(i, :)';
    xact = A*xact + B*uact;
    z = [z; uact; xact];
end
residual = Phi*x0 + Psi*z;                  % find the residual
disp(max(abs(residual)));                   % display the residual

%% Functions
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
function M = otimes(A1, A2)
    M = kron(A1, A2);                       % we can use directly this command as it is programmed already
    return
end