function [F1,F2,F3] = compare_gricsv_with_lqr(sigma)
%COMPARE_GRICSV_WITH_LQR Summary of this function goes here
%   Detailed explanation goes here

clc

F1 = 0;
F2 = 0;
F3 = 0;

A = sigma.A;
B = sigma.B;
%
Q = sigma.Q;
L = sigma.L;
R = sigma.R;

if sigma.discr == false
    sys = ss(A,B,eye(size(A)),0);
else
    sys = ss(A,B,eye(size(A)),0,1);
end

ME1 = [];
try
    [F2,S,E] = lqr(sys,Q,R,L);
    % eig(A-B*K)
catch ME1
    disp(ME1);
end

ME2 = [];
try
    [F1,X,info]=gricsv(sigma);
    disp(info);
    % eig(A+B*F)
catch ME2
    disp(ME2);
end

ME3 = [];
try
    [M,N] = create_hamiltonian_pencil(sigma);
    if sigma.discr == false
        [Xc,ev_c] = gcare(N,M,size(sigma.A,1));
        F3 = -R\(B'*Xc+L');
    else
        [Xd,ev_d] = gdare(N,M,size(sigma.A,1));
        Fd = -(R+B'*Xd*B)\(B'*Xd*A+L');
        F3 = Fd;
    end
catch ME3
    disp(ME3);
end

if sigma.discr == false
    disp('-- GCTARS --');
    %
    if isempty(ME2) && info.has_solution
        [res_gars,res_ars,res_are] = validate_gars(sigma,F1,X,info.V,info.r);
        norms_gricsv = [norm(res_gars),norm(res_ars),norm(res_are)];
        disp('>> Norms of GRICSV experiment');
        disp(norms_gricsv);
        is_stable_gricsv = sum(real(eig(A+B*F1)) < 0) == length(A);
        disp(['Stability of GRICSV solution: ',num2str(is_stable_gricsv)]);
    end
    disp(' ');
    %
    if isempty(ME1)
        [res_gars,res_ars,res_are] = validate_gars(sigma,-F2,S);
        norms_lqr = [norm(res_gars),norm(res_ars),norm(res_are)];
        disp('>> Norms of LQR experiment');
        disp(norms_lqr);
        %        
        is_stable_lqr = sum(real(eig(A-B*F2)) < 0) == length(A);
        disp(['Stability of LQR solution: ',num2str(is_stable_lqr)]);
    end
    disp(' ');
    %
    if isempty(ME3)
        disp('>> Norms of GCARE experiment');
        [s1,s2,s3] = validate_gars(sigma,F3,Xc);
        nv = [norm(s1),norm(s2),norm(s3)];
        disp(nv);
        is_stable = all((sum(real(ev_c) < 0)) == length(A));
        disp(['Stability of GCARE solution: ',num2str(is_stable)]);
    end
else
    % GDTARS
    disp('-- GDTARS --');
    %
    if isempty(ME2) && info.has_solution
        [res_gars,res_ars,res_are] = validate_gars(sigma,F1,X,info.V,info.r);
        norms_gricsv = [norm(res_gars),norm(res_ars),norm(res_are)];
        disp('>> Norms of GRICSV experiment');
        disp(norms_gricsv);
        is_stable_gricsv = sum(abs(eig(A+B*F1)) < 1) == length(A);
        disp(['Stability of GRICSV solution: ',num2str(is_stable_gricsv)]);
    end
    disp(' ');
    %
    if isempty(ME1)
        [res_gars,res_ars,res_are] = validate_gars(sigma,-F2,S);
        norms_lqr = [norm(res_gars),norm(res_ars),norm(res_are)];
        disp('>> Norms of LQR experiment');
        disp(norms_lqr);
        %
        is_stable_lqr = sum((abs(eig(A-B*F2)) < 1)) == length(A);
        disp(['Stability of LQR solution: ',num2str(is_stable_lqr)]);
    end
    disp(' ');
    %
    if isempty(ME3)
        disp('>> Norms of GDARE experiment');
        [s1,s2,s3] = validate_gars(sigma,Fd,Xd);
        nv = [norm(s1),norm(s2),norm(s3)];
        disp(nv);
        is_stable = sum(abs(ev_d) < 1) == length(A);
        disp(['Stability of GDARE solution: ',num2str(is_stable)]);
    end
end

end

