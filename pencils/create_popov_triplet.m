function sigma = create_popov_triplet(A,B,Q,L,R,discr,margin)
%CREATE_POPOV_TRIPLET Declares and validates a Popov tripet sigma=(A,B,P),
% where P = [Q,L;L',R]:
%
%  * For a continuous-time system and performance cost function
%     dx(t)/dt = Ax(t) + Bu(t), with J = Integral{x'Qx + u'Ru + 2*x'Lu}dt
%
%  * For a discrete-time system and performance cost function
%     x[k+1] = Ax[k] + B[k], with J = Sum{x'Qx + u'Ru + 2*x'Lu}.
%
% sigma = create_popov_triplet(A,B,Q,L,R)
%
% sigma = create_popov_triplet(A,B,Q,L,R,discr), where discr==true
%  specifies a discrete-time system, otherwise a continuous-time system
%
% sigma = create_popov_triplet(A,B,Q,L,R,discr,margin)
%

compatible_triplet = ...
    (size(A,1) == size(A,2)) && ...
    (size(B,1) == size(A,1)) && ...
    (size(B,2) == size(R,1)) && ...
    (size(Q,1) == size(Q,2)) && ...
    (size(R,1) == size(R,2)) && ...
    (size(L,1) == size(Q,1)) && ...
    (size(L,2) == size(R,1));
if ~compatible_triplet
    error('Incompatible performance index dimensions.');
end

% if  ~issymmetric(Q)
%     warning('Matrix Q must be symmetric.');
%     Q = Q';
% end
% 
% if  ~issymmetric(R)
%     warning('Matrix R must be symmetric.');
%     R = R';
% end

if nargin == 5
    discr = false;
    margin = -0.2;
elseif nargin == 6
    if discr == false
        margin = -0.2;
    else
        margin = 0.95;
    end
elseif nargin == 7
    
else
    error('Unexpected input arguments.');
end

sigma = struct('A',A,'B',B,'Q',Q,'L',L,'R',R,...
    'discr',discr,'margin',margin);
    
end