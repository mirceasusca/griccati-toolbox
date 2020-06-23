function stab = small_gain_th(G1,G2,gamma)
%SMALL_GAIN_TH Given two systems with Ai stable and the transfer matrices
% G1 and G2 of dimensions p x m and m x p, respectively, the function
% checks if the closed-loop feedback system of G1 and G2 is (internally)
% stable using the small gain theorem.
%
%  The feedback connection is well-posed if and only if the matrix
%   (I  D1) is nonsingular.
%   (D2 I )
%
% stab = small_gain_th(G1,G2), default gamma=1.
% stab = small_gain_th(G1,G2,gamma)
%

if nargin == 2
    gamma = 1;
end

if gamma < 0
    error('Gamma must be positive.');
end

if ~is_stable_ss(G1.a,~(G1.Ts == 0)) || ~is_stable_ss(G2.a,~(G1.Ts == 0))
    error('Input systems G1 and G2 must be stable.');
end

p1 = size(G1.c,1);
m1 = size(G1.b,2);
%
m2 = size(G2.c,1);
p2 = size(G2.b,2);
if p1 ~= p2 || m1 ~= m2
    error('Incompatible system dimensions for feedback connection.');
end

m = m1; p = p1;
Sc = eye(m) - G2.d*G1.d;
Sct = eye(p) - G1.d*G2.d;
if rcond(Sc) < eps || rcond(Sct) < eps
    error('The feedback connection is not well-posed.');
end

stab = false;
if norm(G1,inf) < 1/gamma && norm(G2,inf) <= gamma
    stab = true;
end

CL = feedback(G1,G2);
assert(stab == is_stable_ss(CL.a,~(CL.Ts == 0)));

end