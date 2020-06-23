function [V,W,J,F,X] = kspys_slv(sigma,X)
%KSPYS_SLV Solves the Kalman-Szego-Popov-Yakubovich system in discrete-time
% case.
%
% [V,W,J,F] = kspys_slv(sigma) 
% [V,W,J,F] = kspys_slv(sigma,X)
%

assert(sigma.discr == true,...
    'Function must be called only for discrete-time systems.');

if ~ishermitian(sigma.R)
    warning('R must be hermitian.');
    sigma.R = (sigma.R + sigma.R')/2;
end

A = sigma.A;
B = sigma.B;
Q = sigma.Q;
L = sigma.L;
R = sigma.R;

if nargin == 1
    [~,X,info] = gricsv(sigma);
    has_solution = info.has_solution;
else
    has_solution = true;
end

if has_solution

    Rd = R+B'*X*B;
    Rd = (Rd+Rd')/2;
    [Vs,Ds] = schur(Rd,'real');
    Ds = diag(diag(Ds));
    [Vd,D]=ordschur(Vs,Ds,'lhp');
    J = sign(D);
    
    tol = 1e-12;
    m1 = sum(diag(D) < -tol);
    m2 = sum(diag(D) > tol);
    mz = sum(abs(diag(D)) <= tol);
    assert(mz == 0);
    %
    R11 = Rd(1:m1,1:m1);
    R12 = Rd(1:m1,m1+1:end);
    R22 = Rd(m1+1:end,m1+1:end);
    R11x = R11 - R12/R22*R12';
    %
    V = [sqrtm(-R11x), zeros(size(R12)); R22^(-1/2)*R12', sqrtm(R22)];
    % for debug
    % R+B'*X*B - V'*J*V == 0
    
    W = inv(J)'*inv(V)'*(B'*X*A+L');
    % for debug
    % L+A'*X*B-W'*J*V == 0
    % Q+A'*X*A-X - W'*J*W == 0
    
    F = -V\W;
    % for debug
    % F == -V\W
end

