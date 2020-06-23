function isbounded = bounded_real_lemma(A,B,C,D,gamma,discr)
%BOUNDED_REAL_LEMMA Checks if the continuous-time or discrete-time 
% state-space system (A,B,C,D) and contractness coefficient gamma satisfy 
% the generalized bounded real lemma.
% If attaches the gamma-contractiveness Popov triplet associated with
%      [ A | B ]
%  G = [--- ---], sigma = (A,B;C'*C,C'*D,-gamma^2*I+D'*D,discr).
%      [ C | D ]
% 
% Based on Lemma 7.1.1, the following statements are equivalent:
% I.
%  1) A is stable, ||G||inf < gamma.
%  2) The CTARE(sigma) has a stabilizing solution X >= 0.
% II.
%  1) A is dichotomic, (A,B) stabilizable, ||G|| < gamma.
%  2) The CTARE(sigma) has a stabilizing solution X.
%
% isbounded = bounded_real_lemma(A,B,C,D,gamma)
% isbounded = bounded_real_lemma(A,B,C,D,gamma,discr)
%
% TODO: Corollary 7.1.2: strict bounded real lemma.
%

assert(nargin >= 5,'Not enough input arguments.');
if nargin == 5
    discr = false;
end

if gamma < 0
    error('gamma must be a strictly positive real number.')
end

I = eye(size(B,2));
sigma = create_popov_triplet(A,B,C'*C,C'*D,-gamma^2*I+D'*D,discr);

% for debug
% P = [sigma.Q,sigma.L;sigma.L',sigma.R];
% eig(P)

isbounded = false;

if is_stable_ss(A,discr) % A is stable
    % Part I of the Lemma
    
    [F,X,info] = gricsv(sigma);
    
    if info.has_solution
        semi = true;
        if ispositive(X,semi) && is_stable_ss(A+B*F,discr)
            isbounded = true;
        end
    end

elseif isdichotomic(A,discr)
    % Part II of the Lemma
    
    if isstabilizable(A,B,C,D,discr)
        [F,X,info] = gricsv(sigma);
        if info.has_solution == true
            if is_stable_ss(A+B*F,discr)
                isbounded = true;
            end
        end
    end
    

end

end