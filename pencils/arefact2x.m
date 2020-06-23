function [X,Report] = arefact2x(X1,X2,D)
% Computes Riccati solution X = D*(X2/X1)*D from X1,X2,D.

%   Author(s): Pascal Gahinet
%              edited by Mircea Susca
%   Copyright 1986-2006 The MathWorks, Inc.

% Solve X * X1 = X2
[l,u,p] = lu(X1,'vector');
CondX1 = rcond(u);

if CondX1 > eps
    % Solve for X based on LU decomposition
    X(:,p) = (X2/u)/l;
    
    % Symmetrize
    X = (X+X')/2;
    
    % Factor in scaling D (X -> DXD)
    X = lrscale(X,D,D);
    
    Report = 1;
else
    % X1 is singular
    X = [];  Report = -1;
end

end
