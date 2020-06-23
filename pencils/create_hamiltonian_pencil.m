function [M,N] = create_hamiltonian_pencil(sigma)
%CREATE_HAMILTONIAN_PENCIL Creates the Hamiltonian pencil from the 
% input Popov triplet "sigma".
%
% [M,N] = create_hamiltonian_pencil(sigma)
%
%  * sigma.discr == false, then the system is continuous and declares the
%   extended Hamiltonian pencil (EHP):
%                          (I  0  0)   (A   0   B)
%    lambda*M - N = lambda*(0  I  0) - (Q  -A'  L)
%                          (0  0  0)   (L'  B'  R)
%
%  * sigma.discr == true, then the system is discrete and declares the
%   extended symplectic pencil (ESP):
%                          (I   0   0)   (A   0  B)
%    lambda*M - N = lambda*(0  -A'  0) - (Q  -I  L)
%                          (0  -B'  0)   (L'  0  R)
%

A = sigma.A;
B = sigma.B;
%
n = size(A,1);
m = size(B,2);
%
Q = sigma.Q;
L = sigma.L;
R = sigma.R;
%
M = zeros(2*n+m);
%
if sigma.discr == false
    M(1:2*n,1:2*n) = eye(2*n);
    %
    N = [A, zeros(size(A)), B;
        -Q, -A', -L;
        L', B', R];
else
    M(1:n,1:n) = eye(n);
    M(n+1:2*n,n+1:2*n) = -A';
    M(2*n+1:end,n+1:2*n) = -B';
    %
    N = [A, zeros(size(A)), B;
        Q, -eye(size(A)), L;
        L', zeros(size(B')), R];
end

end