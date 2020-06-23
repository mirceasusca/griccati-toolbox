function [A,B1,B2,C1,C2,D11,D12,D21,D22,p1,m1] = partition_generalized_plant(P,p2,m2)
%PARTITION_GENERALIZED_PLANT Given a generalized plant P in state-space
% form with inputs of size m=m1+m2 and outputs of size p=p1+p2, the 
% function extracts all subblocks as in:
%
% dx/dt =  A*x +  B1*u1 + B2*u2
%    y1 = C1*x + D11*u1 + D12*u2
%    y2 = C2*x + D21*u1 + D22*u2
%
% [A,B1,B2,C1,C2,D11,D12,D21,D22,p1,m1] = partition_generalized_plant(P,p2,m2)
%

A = P.a; B = P.b;
C = P.c; D = P.d;

% check sizes
nx = max(size(A));
[i,m] = size(B);
[p,i] = size(C);
if m <= m2
  disp('control input dimension incorrect')
  return
end
if p <= p2
  disp('measurement output dimension incorrect')
  return
end
p1 = p - p2;
m1 = m - m2;

B1 = B(:,1:m1);
B2 = B(:,m1+1:end);
%
C1 = C(1:p1,:);
C2 = C(p1+1:end,:);
%
D11 = D(1:p1,1:m1);
D12 = D(1:p1,m1+1:end);
D21 = D(p1+1:end,1:m1);
D22 = D(p1+1:end,m1+1:end);

end

