function [V,info] = mpdefsub(M,N,discr,margin,tol,stable)
%MPDEFSUB Computes maximal Cg proper deflating subspace of a matrix pencil.
% Implements Algorithm 2 -- pg. 209.
%
% [V,info] = mpdefsub(M,N)
% [V,info] = mpdefsub(M,N,discr)
% [V,info] = mpdefsub(M,N,discr,margin)
% [V,info] = mpdefsub(M,N,discr,margin,tol)
% [V,info] = mpdefsub(M,N,discr,margin,stable)
%
%   Detailed explanation goes here
% M*V*S = N*S, with S having the desired spectrum (stable/antistable).
%
% See also GSFSTAB
%

info = struct('kron_info',struct(),'ng',-1,'S',0,'nr',0);

[m,n] = size(N);

% for easier interfacing with gklf
A = -N;
E = -M;

% parse input arguments
if nargin == 2
    discr = false;
    margin = -0.2;
    tol = m*n*eps(max(1,norm(A,1)));
    stable = true;
elseif nargin == 3
    if discr == false
        margin = -0.2;
    else
        margin = 0.95;
    end
    tol = m*n*eps(max(1,norm(A,1)));
    stable = true;
elseif nargin == 4
    tol = m*n*eps(max(1,norm(A,1)));
    stable = true;
elseif nargin == 5
    stable = true;
end

% Steps 1) and 2)
[At,Et,kron_info,Q,Z] = gklf(A,E,tol,'reverse'); % (6.59)
info.kron_info = kron_info;

ninf = sum(kron_info.minf);
nf = kron_info.mf;

% Step 3)
[Ar,Er,rkr,rkc] = extract_finite_part(At,Et,kron_info,'right');
info.nr = rkr;
% extract regular part and apply the real Schur form on it
rridx = rkr+1; % starting index for the regular part (row) - right K. str.
rcidx = rkc+1; % starting index for the regular part (column)

[~,~,lkr,lkc] = extract_finite_part(At,Et,kron_info,'left');
% lridx = m-lkr+1; % starting index for the regular part (row) - left K. str
% lcidx = n-lkc+1; % starting index for the regular part (column)

[Ar,Er,Q1,Z1] = qz(Ar,Er,'real');
ev = ordeig(Ar,Er);
% if discr % TODO, for use for numerical refinements
%     select = abs(ev) > 1+tol;
% else
%     select = abs(ev) < 1-tol;
% end
if discr
    select = abs(ev) < 1-tol;
else
    select = (real(ev) < -tol) | (abs(ev)*tol > 1);
end
[Ar,Er,Q1,Z1] = ordqz(Ar,Er,Q1,Z1,select);

ng = sum(select == true);
info.ng = ng; % n1 is the same as ng (6.60)

At(rridx:rridx+nf-1,rcidx:rcidx+nf-1) = Ar; % (6.60)
Et(rridx:rridx+nf-1,rcidx:rcidx+nf-1) = Er;
%
Qa = eye(size(Q));
Qa(rkr+1:rkr+nf,rkr+1:rkr+nf) = Q1';
%
Za = eye(size(Z));
Za(n-lkc-nf-ninf+1:n-lkc-ninf,n-lkc-nf-ninf+1:n-lkc-ninf) = Z1;
%
Q = Q*Qa;
Z = Z*Za;
% for debug
% Q'*A*Z - At == 0
% E - Q*Et*Z' == 0

% Step 4)
Aeps = At(1:rkr,1:rkc);
Eeps = Et(1:rkr,1:rkc);

nz_cols = 0;
for i = 1:size(Eeps,2)
    if all(Eeps(:,i) == 0)
        nz_cols = nz_cols + 1;
    end
end

Aeps1 = Aeps(:,1:nz_cols);
Aeps2 = Aeps(:,nz_cols+1:end);
%
% Eeps1 = Eeps(:,1:nz_cols); % all zeros, not needed
Eeps2 = Eeps(:,nz_cols+1:end);

opt = struct('stable',stable,'discr',discr);
F = gsfstab(Aeps2,Eeps2,Aeps1,[],margin,opt);
% for debug
% eig(Aeps2,Eeps2)
% eig(Aeps2+Aeps1*F,Eeps2) % as returned from gsfstab

assert(size(Aeps2,1) == size(Aeps2,2));
assert(isempty(F) == (rkr == 0));

if ~isempty(F)
    Z4 = eye(size(Z));
    Z4(1:size(F,1),...
       size(Aeps1,2)+1:size(Aeps1,2)+size(F,2)) = F;
    %
    Z = Z*Z4;
    At = At*Z4;
    Et = Et*Z4;
    
    % apply column permutations -- switch 1
    At(:,1:size(Aeps,2)) = At(:,[size(Aeps1,2)+1:size(Aeps,2),1:size(Aeps1,2)]);
    Et(:,1:size(Aeps,2)) = Et(:,[size(Aeps1,2)+1:size(Aeps,2),1:size(Aeps1,2)]);
    Z5 = eye(size(Z));
    Z5(:,1:size(Aeps,2)) = Z5(:,[size(Aeps1,2)+1:size(Aeps,2),1:size(Aeps1,2)]);
    Z = Z*Z5;

    % apply column permutations -- switch 2
    At(:,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng) = At(:,[size(Aeps2,2)+size(Aeps1,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)]);
    Et(:,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng) = Et(:,[size(Aeps2,2)+size(Aeps1,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)]);
    Z6 = eye(size(Z));
    Z6(:,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng) = Z6(:,[size(Aeps2,2)+size(Aeps1,2)+1:size(Aeps2,2)+size(Aeps1,2)+ng,size(Aeps2,2)+1:size(Aeps2,2)+size(Aeps1,2)]);
    Z = Z*Z6;
end

ns = rkr + ng;
V = Z(:,1:ns);

A1 = Ar(1:ng,1:ng);
E1 = Er(1:ng,1:ng);
if ~isempty(F) % singular pencil
    X2 = At(1:size(Aeps2,1),size(Aeps2,2)+1:size(Aeps2,2)+ng);
    X1 = Et(1:size(Aeps2,1),size(Aeps2,2)+1:size(Aeps2,2)+ng);
    %
    S11 = Eeps2\(Aeps1*F+Aeps2);
    S12 = Eeps2\(X2-X1/E1*A1);
    S22 = E1\A1;
    %
    S = zeros(length(S11)+length(S22));
    S(1:size(S11,1),1:size(S11,2)) = S11;
    S(1:size(S11,1),size(S11,2)+1:end) = S12;
    S(size(S11,1)+1:end,size(S11,2)+1:end) = S22;
else % regular pencil
    S = E1\A1;
end
info.S = S;

end

function [Ar,Er,krd,krc] = extract_finite_part(At,Et,info,structure)
%EXTRACT_REGULAR_PART
% krd - kronecker row dimension
% krc - kronecker column dimension

if nargin < 4
    structure = 'right';
end

if strcmp(structure,'right')
    ridx = 1;
    cidx = 1;
    if ~isempty(info.mr)
        for i = 1:length(info.mr)
            ridx = ridx + (i-1)*(info.nr(i)-info.mr(i));
            cidx = cidx + i*(info.nr(i)-info.mr(i));
        end
    end

    nf = info.mf;

    Ar = At(ridx:ridx+nf-1,cidx:cidx+nf-1);
    Er = Et(ridx:ridx+nf-1,cidx:cidx+nf-1);

    krd = ridx-1;
    krc = cidx-1;
elseif strcmp(structure,'left')
    ridx = 1;
    cidx = 1;
    if ~isempty(info.ml)
        ll = length(info.ml);
        for i = ll:-1:1
            cidx = cidx + (ll-i)*(info.ml(i)-info.nl(i));
            ridx = ridx + (ll-(i-1))*(info.ml(i)-info.nl(i));
        end
    end

    Ar = 0;
    Er = 0; % don't care for left call

    krd = ridx-1;
    krc = cidx-1;
else
    error('Unknown requested structure.');
end

end
