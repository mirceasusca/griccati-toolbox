function hassol = has_pos_sd_stab_sol(A,B,C,D,discr)
%HAS_POS_SD_STAB_SOL Checks if the ARE (CTARE/DTARE) of the positiveness 
% Popov triplet associated with the system (A,B,C,D) has a positive 
% semidefinite stabilizing solution X.
%
% Uses Proposition 7.2.1, pg. 219.
%
% hassol = has_pos_sd_stab_sol(A,B,C,D)
% hassol = has_pos_sd_stab_sol(A,B,C,D,discr)
%

if nargin == 4
    discr = false;
end

sigma = create_popov_triplet(A,B,C'*C,C'*D,D'*D,discr);

% preconditions
% a) (A,B) is stabilizable
c1 = isstabilizable(A,B,discr);

if discr == false
    % b) D is of full column rank
    c2 = rank(D) == size(D,2);
else
    % no need to check this condition for the discrete-time case
    c2 = true;
end

% c) [jwI-A, -B] is of full column rank for all w in R, for a CT system.
%    [    C,  D]
%
%  [exp(j*theta)I-A, -B] is of full column rank for all theta in [0,2*pi),
%  [              C,  D] for a DT system.
M = [eye(size(A)), zeros(size(B)); zeros(size(C)), zeros(size(D))];
N = [A,B;-C,-D];
tol = size(N,1)*size(N,2)*eps(max(1,norm(N,1)));
% NI(2) contains the normal rank of the system pencil
[zf,ni] = sl_gzero(N,M,[],[],[],tol);
c3 = ni(2) == size(M,2);

hassol = c1 && c2 && c3;

% for debugging
if hassol
    [F,X,info] = gricsv(sigma);
    if info.has_solution
        semi = true;
        assert((info.has_solution == true) && ispositive(X,semi));
    else
        warning(['The function returned that the ARE has a solution',...
            ' and GRICSV failed to find it.']);
    end
else
    [F,X,info] = gricsv(sigma);
    assert(info.has_solution == false);
end

end
