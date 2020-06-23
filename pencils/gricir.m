function X1 = gricir(Sigma,X)
%GRICIR Refines the solution X to the Algebraic Riccati Equations
% for the continuous or discrete-time Popov triplet Sigma using a Newton
% iteration for the attached Lyapunov equation.
%
% Bibliography for the continuous-time case with L = 0:
% Numerical Methods for Linear Control Systems, Biswa Nath Datta, 
% Ch. 13 -- NUMERICAL SOLUTIONS AND CONDITIONING OF ALGEBRAIC
% RICCATI EQUATIONS.
%
% X_new = gricir(Sigma,X)
%

A = Sigma.A;
B = Sigma.B;
Q = Sigma.Q;
L = Sigma.L;
R = Sigma.R;

if Sigma.discr == false
    RCC = @(X) A'*X+X*A-(X*B+L)/R*(X*B+L)'+Q;
    S = R\(X*B+L)';
    Ak = A-B*S;
    dXk = lyap(Ak',RCC(X));
    X1 = X + dXk;
else
    RCD = @(X) A'*X*A-X+Q-(L+A'*X*B)/(R+B'*X*B)*(L'+B'*X*A);
    K = (R+B'*X*B)\(B'*X*A+L');
    Ak = A-B*K;
    dXk = dlyap(Ak',RCD(X));
    X1 = X + dXk;
end
