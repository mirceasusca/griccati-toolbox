function [F,X,info] = gricsv(Sigma,bal)
%GRICSV Solves the generalized algebraic Riccati system (GCTARS/GDTARS)
% described by the input Popov triplet "Sigma".
% Implements Algorithm 3 -- pg. 212.
%
% [F,X,info] = gricsv(sigma)
% [F,X,info] = gricsv(sigma,'balance')
% [F,X,info] = gricsv(sigma,'nobalance')
%
%   Detailed explanation goes here
%
% See also MPDEFSUB
%

if nargin == 1
    bal = 'balance';
end

info = struct('has_solution',false,'r',0,...
    'system','GCTARS','regular',false,'V',0,'S',0);

F = Inf; X = Inf;

if Sigma.discr == true
    info.system = 'GDTARS';
end

n = size(Sigma.A,1); % system order
m = size(Sigma.B,2); % number of inputs

[M,N] = create_hamiltonian_pencil(Sigma);
% for debug
% M0 = M; N0 = N;
% [N1,M1,sx,perm] = arescale(N,M,n);
% D = eye(n);
% perm = 1;

if strcmp(bal,'balance')
    [N,M,sx,perm] = arescale(N,M,n);
    D = diag(sx);
    if ~all(sx == 1)
        disp('Applied balancing.');
    end
elseif strcmp(bal,'nobalance')
    D = eye(n);
    sx = ones(1,n);
    perm = 1:n;
else
    error('Unaccepted value for parameter bal.');
end

% compute tolerance
tol = size(N,1)*size(N,2)*eps(max(1,norm(N,1)));

stable = true;

% numerical refinements
% revert matrices to avoid reordering of eigenvalues
% M0 = M; % also, change selection for reordering in MPDEFSUB
% N0 = N;
% if sigma.discr
%     M = N0;
%     N = M0;
%     sigma.margin = 1.05;
%     stable = false;
% else
%     M = N0 + M0;
%     N = N0 - M0;
%     sigma.margin = 0.95;
%     stable = true;
% end

% Step 1
[V,mpdef_info] = mpdefsub(M,N,Sigma.discr,Sigma.margin,tol,stable);
V(abs(V) < tol) = 0;
[rnk,nu_r,nu_l] = nrankpi(mpdef_info.kron_info,size(M,1),size(M,2));
info.V = V;
info.S = mpdef_info.S;

nr = mpdef_info.nr;
ng = mpdef_info.ng;

% if nu_r > 0
if rnk ~= size(M,1)
    info.regular = false;
    warning('Pencil is singular.');
else
    info.regular = true;
end

% Step 2)
if nr + ng < n
    warning(['There exists no solution of maximal dimension to the ',info.system,'.']);
    return;
end

% Step 3) already done in mpdefsub
% if nr == 0
%     % check V for the regular case computed using Steps 1-3 of algorithm 2
%     disp('No right Kronecker structure.');
% else
%     disp('Has right Kronecker structure.');
%     % check V computed using Step 4 of algorithm 2
% end

% Step 4) partition V as in (6.6) or (6.31)
V1 = V(1:n,:);
V2 = V(n+1:2*n,:);
V3 = V(2*n+1:end,:);

if rank(V1) < n
    warning(['There exists no solution of maximal dimension to the ',info.system,'.']);
    return;
end

% Step 5)
r = nr + ng;
info.r = r;

assert(r <= n,'r > n, check computations and tolerances.');
% if ~all(perm == (1:length(perm))')
%     warning('Permutation matrix not unitary.');
% end
    
disp(['Condition number of V1: ',num2str(cond(V1))]);
if r == n
    % CTARS/DTARS solution: (X,F)
    X = arefact2x(V1,V2,sx);
    X = gricir(Sigma,X);

    % method 1) - the best solution, as it has the lowest residual R(X)->0
    if Sigma.discr == false
        F = -Sigma.R\(Sigma.B'*X+Sigma.L');
    else
        F = -(Sigma.R+Sigma.B'*X*Sigma.B)\(Sigma.B'*X*Sigma.A+Sigma.L');
    end
    
    % method 2)
    % [V,W,J,F] = kpys(Sigma,X);
    
    % method 3)
    % F = V3/V1*D;
elseif r < n
    % GCTARS/GDTARS solution: (X,F,r,V1)
    V1pinv = pinv(V1);
    X = D*(V2*V1pinv)*D;
    X = (X+X')/2; % symmetrize solution
    %
    X = gricir(Sigma,X);

    % method 1)
    if Sigma.discr == false 
        F = -Sigma.R\(Sigma.B'*X+Sigma.L');
    else
        F = -(Sigma.R+Sigma.B'*X*Sigma.B)\(Sigma.B'*X*Sigma.A+Sigma.L');
    end
    
    % method 2)
    % [V,W,J,F] = kpys(Sigma,X);
    
    % method 3)
    % F = V3*V1pinv*D;
else
    warning(['This case should not be achievable for the ',info.system,'.']);
end
info.has_solution = true;


end
