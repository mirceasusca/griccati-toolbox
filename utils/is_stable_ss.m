function isst = is_stable_ss(A,discr)
%IS_STABLE_SS Checks if the state matrix A is stable in continuous-time or
% discrete-time, based on the argument discr.
%
% isst = is_stable_ss(A)
% isst = is_stable_ss(A,discr)
%

if size(A,1) ~= size(A,2)
    error('Matrix must be square.');
end

if nargin == 1
    discr = false;
end

if discr == false
    isst = isnegative(A);
else
    ev = eig(A);
    tol = size(A,1)*size(A,1)*eps(max(1,norm(A,1)));
    isst = all(abs(ev) < 1 - tol);
end

end