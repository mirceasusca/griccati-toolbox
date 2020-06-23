function [K,CL,GAM,INFO] = hinfgsyn(P,p2,m2,gam,bal)
%HINFGSYN Solves the suboptimal H-infinity control problem.
%  Searches for a proper controller K for which the closed-loop system 
% T(y1,u1) = LFT(T,K) is internally stable, i.e. lambda(AR) in C- and 
% has the H-infinity norm lower than the prescribed input value "gam". 
% Otherwise, returns GAM=Inf and K and CL are set to [].
%  The function computes the controller K using the two-Riccati equations
% method, solved using deflating subspaces.
%
% The plant T is partitioned as:
%                 [A  |  B1  B2]
% T = [T11 T12] = [------------], with
%     [T21 T22]   [C1 | D11 D12]
%                 [C2 | D21 D22]
% u1 - external inputs (disturbance, sensor noise, reference signals)
% u2 - controlled inputs
% y1 - regulated outputs (error signals)
% y2 - measured outputs
%
% [K,CL,GAMMA,INFO] = HINFGSYN(PLANT,NMEAS,NCONT,GAMTRY)
% [K,CL,GAMMA,INFO] = HINFGSYN(PLANT,NMEAS,NCONT,GAMTRY,'balance')
% [K,CL,GAMMA,INFO] = HINFGSYN(PLANT,NMEAS,NCONT,GAMTRY,'nobalance')
%
% The following hypotheses hold:
% (h1) (A,B2) stabilizable and (C2,A) detectable;
% (h2) Regularity
%   (r1) rank(D12) = m2, rank [jwI-A,  -B2] = n + m2, for all w in R
%                             [  -C1, -D12];
% 
%   (r2) rank(D21) = m2, rank [jwI-A,  -B1] = n + m2, for all w in R
%                             [  -C2, -D21].
%

if nargin == 4
    bal = 'balance';
end
if ~(strcmp(bal,'balance') || strcmp(bal,'nobalance'))
    error('Unaccepted value for parameter bal. Check help for usage.');
end

K = [];
CL = [];
GAM = Inf;
INFO = struct('Xinfo',[],'Yinfo',[],'Zinfo',[]);

P = scale_gen_plant_gamma(P,p2,m2,gam);
[A,B1,B2,C1,C2,D11,D12,D21,D22,p1,m1] = partition_generalized_plant(P,p2,m2);
discr = (P.Ts ~= 0);

% it follows that p1>=m2, i.e. there are more regulated outputs than
% controlled inputs, and T12 has no zeros on C0bar
c1 = has_pos_sd_stab_sol(A,B2,C1,D12,discr);

% it follows that m1>=p2, i.e. there are more external inputs than
% measured outputs, and T21 has no zeros on C0bar
c2 = has_pos_sd_stab_sol(A',C2',B1',D21',discr);

if ~(c1 && c2)
    warning(['The Hinf control problem is singular and',...
        ' may not have a solution. Try regularizing the problem.']);
end

% the solution
Qc = C1'*C1;
Lc = C1'*[D11,D12];
Rc = [D11';D12']*[D11,D12] - ...
    blkdiag(eye(size(D11,2)),zeros(size(D12,2)));
sigma_c = create_popov_triplet(A,[B1,B2],Qc,Lc,Rc,discr);
%
Qo = B1*B1';
Lo = B1*[D11',D21'];
Ro = [D11;D21]*[D11',D21'] - ...
    blkdiag(eye(size(D11,1)),zeros(size(D21*D21')));
sigma_o = create_popov_triplet(A',[C1',C2'],Qo,Lo,Ro,discr);

[Fc,X,info_c] = gricsv(sigma_c,bal);
[Fo,Y,info_o] = gricsv(sigma_o,bal);
INFO.Xinfo = info_c;
INFO.Yinfo = info_o;

% checking (I)
c1 = info_c.has_solution && ispositive(X,true);
c2 = info_o.has_solution && ispositive(Y,true);
if ~(c1 && c2)
    warning('Solutions X and Y are not positive semidefinite.');
    return;
end

rho = max(abs(eig(X*Y)));
c3 = rho < 1;
if c3 == false
    warning('Spectral radius of X*Y is not strictly subunitary.');
    return;
end

% checking (II) % the Riccati equations have positive semidefinite
% solutions iff the KPYS have positive semidefinite solutions (X,V,W)
[Vc,Wc,Jc] = kpys(sigma_c,X);
[Vo,Wo,Jo] = kpys(sigma_o,Y);

% partition KPYS(sigma_c,Jc)
Vc11 = Vc(1:m1,1:m1);
Vc12 = Vc(1:m1,m1+1:end); % only for assertion
Vc21 = Vc(m1+1:end,1:m1);
Vc22 = Vc(m1+1:end,m1+1:end);
assert(all(Vc12(:) < eps),'Vc not lower block triangular.')
% Wc1 = Wc(1:m1,:);
% Wc2 = Wc(m1+1:end,:);
%
% Jx = blkdiag(-eye(m2),eye(p2));
Qx = B1/(Vc11'*Vc11)*B1';
Lx = B1/(Vc11'*Vc11)*[Vc21',D21'];
Rx = [Vc21;D21]/(Vc11'*Vc11)*[Vc21',D21'] - blkdiag(eye(m2),zeros(size(D21*D21')));
%
F1 = Fc(1:m1,:);
F2 = Fc(m1+1:end,:); % length m2

% build KPYS(sigma_x,Jx)
sigma_x = create_popov_triplet(A'+F1'*B1',[-F2'*Vc22',C2'+F1'*D21'],Qx,Lx,Rx,discr);

[Fx,Z,info_z] = gricsv(sigma_x);
INFO.Zinfo = info_z;
[Vx,Wx,Jx] = kpys(sigma_x); % kpys(sigma_x,Z);
Vx11 = Vx(1:m2,1:m2);
Vx12 = Vx(1:m2,m2+1:end); % only for assertion
% Vx21 = Vx(m2+1:end,1:m2);
% Vx22 = Vx(m2+1:end,m2+1:end);
assert(all(Vx12(:) < eps),'Vx not lower block triangular.')

% Z = Y/(eye(size(X*Y))-X*Y); % more accurate than using the KPYS(sigma)
% assert(norm(Z - Y/(eye(size(X*Y))-X*Y),1) < 1e-3);

C2F1 = C2 + D21*F1;
% Sc = inv(Vc11'*Vc11);
Sx = inv(Vx11'*Vx11);

if discr == false

    Dg11 = -((D12'*D12)\D12')*D11/(Vc11'*Vc11)*(D21'/(D21/(Vc11'*Vc11)*D21'));
    Dg21 = sqrtm(inv(D21/(Vc11'*Vc11)*D21'));
    Cg2 = ((D21/(Vc11'*Vc11)*D21')^(-1/2))*C2F1;
    Bg1 = -(B1/(Vc11'*Vc11)*D21')/(D21/(Vc11'*Vc11)*D21') - B2*Dg11 - Z*Cg2'*Dg21;
    Ag = A + B1*F1 + B2*F2 + Bg1*C2F1;
    Cg1 = -F2 + Dg11*C2F1;
    Dg12 = (D12'*D12)^(-1/2)*Sx^(-1/2);
    Bg2 = -B1/(Vc11'*Vc11)*((D12'*D12)\D12'*D11 + ...
        Dg11*D21)'*sqrtm(D12'*D12)*sqrtm(Sx)...
        - B2*Dg12 - Z*Cg1'*sqrtm(D12'*D12)*sqrtm(Sx);
    Dg22 = zeros(size(Dg21,1),size(Dg12,2));

else

    % Vc22^-1*Vc21 = (D12'*D12+B2'*X*B2)\(D12'*D11+B2'*X*B1);
    % Vc22 = Vc22' = sqrtm(D12'*D12+B2'*X*B2)
    
    Dg11 = -((D12'*D12+B2'*X*B2)\(D12'*D11+B2'*X*B1)/(Vc11'*Vc11)*D21' - F2*Z*C2F1')/(D21/(Vc11'*Vc11)*D21'+C2F1*Z*C2F1');
    Dg21 = (D21/(Vc11'*Vc11)*D21'+C2F1*Z*C2F1')^(-1/2);
    Cg2 = ((D21/(Vc11'*Vc11)*D21'+C2F1*Z*C2F1')^(-1/2))*C2F1;
    Bg1 = -(B1/(Vc11'*Vc11)*D21')/(D21/(Vc11'*Vc11)*D21'+C2F1*Z*C2F1') - B2*Dg11 - ...
        (A+B1*F1)*Z*Cg2'*Dg21;
    Ag = A + B1*F1+B2*F2+Bg1*C2F1;
    Cg1 = -F2 + Dg11*C2F1;
    Dg12 = Vc22\(Sx^(-1/2));
    Bg2 = -B2*Dg12 - B1/(Vc11'*Vc11)*((D12'*D12+B2'*X*B2)\(D12'*D11+B2'*X*B1) + Dg11*D21)'*sqrtm(D12'*D12+B2'*X*B2)'*sqrtm(Sx) - ...
        (A+B1*F1)*Z*Cg1'*sqrtm(D12'*D12+B2'*X*B2)'*sqrtm(Sx);
    Dg22 = zeros(size(Dg21,1),size(Dg12,2));
    
end

% % for debug
% AK = Ag;
% BK = Bg1;
% CK = Cg1;
% DK = Dg11;
% % K = ss(AK,BK,CK,DK,P.Ts); % central controller
% assert(size(AK,1) == size(AK,2))
% assert(size(BK,1) == size(AK,1))
% assert(size(BK,2) == p2);
% assert(size(CK,1) == m2);
% assert(size(CK,2) == size(AK,2))
% assert(all(size(DK) == [m2,p2]))
% %
% S1 = eye(p2)-D22*DK;
% S2 = eye(m2)-DK*D22;
% assert(rank(S1) == length(S1));
% assert(rank(S2) == length(S2));

% Kg = ss(Ag,[Bg1,Bg2],[Cg1;Cg2],[Dg11,Dg12;Dg21,Dg22],P.Ts); % unscaled
Kg = scale_hinf_controller_gamma(Ag,Bg1,Bg2,Cg1,Cg2,...
    Dg11,Dg12,Dg21,Dg22,P.Ts,1/gam);

Q = zeros(m2,p2); % free parameter with ||Q||_inf < gamma
K = lft(Kg,Q);

P = scale_gen_plant_gamma(P,p2,m2,1/gam);

if ~all(D22(:) == 0)
    K = feedback(K,P.d(p1+1:end,m1+1:end));
end

CL = lft(P,K);
INFO.ev = eig(CL.a);
GAM = norm(CL,inf);

end

