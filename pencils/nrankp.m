function [r,nu_r,nu_l] = nrankp(A,E,tol)
%NRANKP Computes normal rank of matrix pencil A-lambda*E.
% nu_r := n - r > 0 then the pencil has right singular structure
% nu_l := m - r > 0 then the pencil has left singular structure
%
% [r,nu_r,nu_l] = nrankp(A,E,tol)
%

[m,n] = size(A);
[m1,n1] = size(E);

if m ~= m1 || n ~= n1
    error('Incompatible matrix pencil.');
end

if nargin < 3
    tol = size(A,1)*size(A,1)*eps(max(1,norm(A,1)));
end

[~,~,info] = gklf(A,E,tol);

[r,nu_r,nu_l] = nrankpi(info,m,n);

end

