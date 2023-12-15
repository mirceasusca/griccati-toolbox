function [Pscaled] = scale_gen_plant_gamma(P,p2,m2,gamma)
%SCALE_GEN_PLANT_GAMMA Conversion of a gamma-DAP with arbitrary gamma > 0 
% to a DAP (with gamma = 1).
%
% Pscaled = scale_gen_plant_gamma(P,p2,m2,gamma)
%

if gamma < 0
    error('Gamma must be strictly greater than zero.');
end

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

B1 = gamma^(-1/2)*B(:,1:m1);
B2 = gamma^(+1/2)*B(:,m1+1:end);
%
C1 = gamma^(-1/2)*C(1:p1,:);
C2 = gamma^(+1/2)*C(p1+1:end,:);
%
D11 = gamma^(-1)*D(1:p1, 1:m1);
D12 = D(1:p1,m1+1:end);
D21 = D(p1+1:end,1:m1);
D22 = gamma*D(p1+1:end,m1+1:end);

Pscaled = ss(A,[B1,B2],[C1;C2],[D11,D12;D21,D22],P.Ts);

end

