function sys1 = inv_sys(sys)
%INV_SYS

A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;

A1 = A - B/D*C;
B1 = B/D;
C1 = -D\C;
D1 = inv(D);

sys1 = ss(A1,B1,C1,D1,sys.Ts);

end