function [r,nu_r,nu_l] = nrankpi(info,m,n)
%NRANKPI Computes normal rank of matrix pencil A-lambda*E given the info 
% stucture returned from a gklf call.
%
% nu_r := n - r > 0 then the pencil has right singular structure
% nu_l := m - r > 0 then the pencil has left singular structure
%
% [r,nu_r,nu_l] = nrankpi(A,E,tol)
%

r = sum(info.minf) + info.mf;
[rkr,rkc] = extract_kronecker_structure(info,'right');
[lkr,lkc] = extract_kronecker_structure(info,'left');
r = r + rkr + lkc;

nu_r = n - r;
nu_l = m - r;

assert(r <= m && r <= n);

end

