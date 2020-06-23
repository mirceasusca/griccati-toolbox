function [r_gars,r_ars,r_are] = validate_gars(sigma,F,X,V,r)
%CTARS Summary of this function goes here
%   Detailed explanation goes here

A = sigma.A;
B = sigma.B;
%
discr = sigma.discr;
%
Q = sigma.Q;
L = sigma.L;
R = sigma.R;

if nargin == 3
    generalized = false;
elseif nargin == 5
    generalized = true;
else
    error('Must receive arguments (sigma,F,X) or (sigma,F,X,V,r).');
end

n = size(A,1);

r_gars = Inf;
if discr == false
    if generalized == true
        r_gars = [A'*X+X*A+Q,L+X*B;B'*X+L',R]*[eye(size(A)); F]*V(1:n,1:r);
    end
    r_ars = [A'*X+X*A+Q, L+X*B; B'*X+L', R]*[eye(size(A)); F];
    r_are = A'*X+X*A + (X*B+L)*F + Q;
else
    if generalized == true
        r_gars = [A'*X*A-X+Q,L+A'*X*B; B'*X*A+L',R+B'*X*B]*...
            [eye(size(A));F]*V(1:n,1:r);
    end
    r_ars = [A'*X*A-X+Q,L+A'*X*B;B'*X*A+L',R+B'*X*B]*[eye(size(A));F];
    r_are = A'*X*A - X + (A'*X*B+L)*F + Q;
end
    
end

