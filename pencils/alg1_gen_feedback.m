function F = alg1_gen_feedback(M,N,B)
% from GRT, pencil = M*lambda - N

[regular,Sr,Tr,Qr,Zr,kstr]=isregular(M,N); % pencil = N - lambda M
% Q'*M*Z-Mr
% Q'*N*Z-Nr
if ~regular
    error('Pencil must be regular.');
end

[S,T,QQ,ZZ]=part_reg_inf_n1_n2(Sr,Tr,Qr,Zr);
Q = Qr*QQ;
Z = ZZ*Zr;
% QQ'*Q'*M*Z*ZZ'-S
% QQ'*Q'*N*Z*ZZ'-T

F = 0;