function [K,CL,GAM,INFO] = h2gsyn(P,p2,m2,bal)
%H2GSYN Solves the optimal H2 control problem.
%  Searches for a (strictly) proper controller K for which the closed-loop
% system T(y1,u1) = LFT(T,K) is internally stable, i.e. lambda(AR) in C- 
% and has minimum H2 norm, i.e. norm(T(y1,u1),2) attains its minimum over 
% the class of systems K = {K: K strictly proper and T(y1,u1) internally stable}.
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
% y1 - regulated outputs (error signals)
% u2 - controlled inputs
% y2 - measured outputs
%
% [K,CL,GAMMA,INFO] = H2GSYN(PLANT,NMEAS,NCONT)
% [K,CL,GAMMA,INFO] = H2GSYN(PLANT,NMEAS,NCONT,'balance')
% [K,CL,GAMMA,INFO] = H2GSYN(PLANT,NMEAS,NCONT,'nobalance')
%
% Assumptions:
% (h1) T11(inf) = D11 = 0;
% (h2) (A,B2) stabilizable and (C2,A) detectable;
% (h3) Regularity
%   (r1) rank(D12) = m2, rank [jwI-A,  -B2] = n + m2, for all w in R
%                             [  -C1, -D12];
% 
%   (r2) rank(D21) = m2, rank [jwI-A,  -B1] = n + m2, for all w in R
%                             [  -C2, -D21].
%

if nargin == 3
    bal = 'balance';
end
if ~(strcmp(bal,'balance') || strcmp(bal,'nobalance'))
    error('Unaccepted value for parameter bal. Check help for usage.');
end

K = [];
CL = [];
GAM = Inf;
INFO = struct('Xinfo',[],'Yinfo',[]);

[A,B1,B2,C1,C2,D11,D12,D21,D22] = partition_generalized_plant(P,p2,m2);
discr = (P.Ts ~= 0);

sigma12 = create_popov_triplet(A,B2,C1'*C1,C1'*D12,D12'*D12,discr);
sigma21 = create_popov_triplet(A',C2',B1*B1',B1*D21',D21*D21',discr);

% test hypotheses
c1 = has_pos_sd_stab_sol(A,B2,C1,D12,discr);
c2 = has_pos_sd_stab_sol(A',C2',B1',D21',discr);

if ~(c1 && c2)
    warning(['The H2 control problem is singular and',...
        ' may not have a solution. Try regularizing the problem.']);
end

% the solution
[Fs,Xs,Xinfo] = gricsv(sigma12,bal);
[Ks,Ys,Yinfo] = gricsv(sigma21,bal);
INFO.Xinfo = Xinfo;
INFO.Yinfo = Yinfo;

if Xinfo.has_solution && Yinfo.has_solution
    % K = ss(A+B2*Fs+Ks'*C2+B2'*D22*C2,-Ks',Fs,0,P.Ts);
    K = ss(A+B2*Fs+Ks'*C2,-Ks',Fs,0,P.Ts);

    if ~all(D22(:) == 0)
        K = feedback(K,D22);
    end

    CL = lft(P,K);
    INFO.ev = eig(CL.a);

    if is_stable_ss(CL.a,discr)
        GAM = norm(CL,2);
        % Lc = gram(CL,'c');
        % Lo = gram(CL,'o');
        % GAM = sqrt(trace(CL.c*Lc*CL.c'));
        % GAM = sqrt(trace(CL.b'*Lo*CL.b));
        % GAM = sqrt(trace(B1'*Xs*B1) + trace(B2'*Xs*Ys*Xs*B2)); % !!
        % GAM = sqrt(trace(B1'*Xs*B1) + trace(C1*Ys*C1'));
    else
        GAM = Inf;
    end
else
    return;

end
