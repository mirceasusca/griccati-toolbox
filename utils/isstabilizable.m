function iss = isstabilizable(A,B,C,D,discr)
%ISSTABILIZABLE Checks stabilizability of the state-space system (A,B,C,D).
% If the argument discr is true, the algorithm checks stabilizability for
% the discrete system (A,B,C,D,Ts=1).
%
% iss = isstabilizable(A,B)
% iss = isstabilizable(A,B,discr)
% iss = isstabilizable(A,B,C,D)
% iss = isstabilizable(A,B,C,D,discr)
%
% See also ISDETECTABLE

if nargin == 4 || nargin == 2
    discr = false;
else
    if nargin == 3
        discr = C;
    end    
end

if nargin == 2 || nargin == 3
    C = eye(size(A));
    D = 0;
end

iss = false;

if discr == false
    G = ss(A,B,C,D);
else
    G = ss(A,B,C,D,1); % Ts = 1 not relevant
end

[GS,GNS] = stabsep(G);
[Abar,Bbar,Cbar,T,K] = ctrbf(GNS.a,GNS.b,GNS.c);

if sum(K) == size(GNS.a,1) % is stabilizable
    iss = true;
end

end

