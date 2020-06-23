function neg = isnegative(X,semi)
%ISNEGATIVE Checks if input square matrix is negative definite.
%
% neg = isnegative(X)
% neg = isnegative(X,semi)
%  semi: if true, checks for negative semidefiniteness (default false)
%
% See also ISPOSITIVE

% check validity for semi
if nargin == 1
    semi = false;
end

neg = ispositive(-X,semi);

end