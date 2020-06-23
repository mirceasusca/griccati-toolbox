% uncontrollable system -- LQR returns an unstabilizing feedback
A = [1,0;0,1];
B = [1;2];
Q = eye(2); 
L = [0;0];
R = 0.01;
sigma = create_popov_triplet(A,B,Q,L,R);

[K,S,E]=lqr(A,B,Q,R,L);
eig(A-B*K)

controllability_rank = rank(ctrb(A,B))

[F,X,info] = gricsv(sigma);
if info.has_solution
    [s1,s2,s3] = validate_gars(sigma,F,X,info.V,info.r);
    eig(A+B*F)
end
