function isd = isdetectable(A,B,C,D,discr)
%ISDETECTABLE Checks detectability of the state-space system (A,B,C,D).
% If the argument discr is true, the algorithm checks detectability for
% the discrete system (A,B,C,D,Ts=1).
%
% isd = isdetectable(A,B)
% isd = isdetectable(A,B,discr)
% isd = isdetectable(A,B,C,D)
% isd = isdetectable(A,B,C,D,discr)
%
% See also ISSTABILIZABLE
%

if nargin == 2 || nargin == 3
    C = eye(size(A));
    D = 0;
end

if nargin == 4 || nargin == 2
    discr = false;
end

isd = isstabilizable(A',C',B',D',discr);

end