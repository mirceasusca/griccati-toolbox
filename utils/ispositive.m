function pos = ispositive(X,semi)
%ISPOSITIVE Checks if input square matrix is positive definite.
%
% pos = ispositive(X)
% pos = ispositive(X,semi)
%  semi: if true, checks for positive semidefiniteness (default false)
%
% See also ISNEGATIVE

% check validity for semi
if nargin == 1
    semi = false;
end

pos = false;
ev = eig(X);

tol = size(X,1)^2*eps(max(1,norm(X,1)));

if semi == true
    if all(real(ev) >= -tol) % uncertainty in the interval [-tol,tol]
        pos = true;
    end
else
    if all(real(ev) > tol)
        pos = true;
    end
end

end
