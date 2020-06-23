function isd = isdichotomic(X,discr)
%ISDICHOTOMIC Checks if square input matrix is dichotomic, i.e. it does not
% have the eigenvalues on the stability limit boundary. For continuous-time
% systems, the boundary represents the imaginary axis, while for
% discrete-time systems, the boundary represents the unit circle.
%
% isd = isdichotomic(X)
% isd = isdichotomic(X,discr), with discr=false as default
%

if size(X,1) ~= size(X,2)
    error('Matrix must be square.');
end

if nargin == 1
    discr = true;
end

tol = size(X,1)^2*eps(max(1,norm(X,1)));
ev = eig(X);

if discr == false
    isd = all(abs(ev) > tol);
else
    isd = all(abs(abs(ev)-1.0) > tol);
end

end

