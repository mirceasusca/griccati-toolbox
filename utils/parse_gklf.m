function parse_gklf(A,E,tol)
%PARSE_GKLF Calls gklf(A,E) for the matrix pencil A - lambda*E and displays
% the blocks of the Kronecker-like structure. TOL is used to determine the
% rank of the individual blocks.
%  
% parse_gklf(A,E)
% parse_gklf(A,E,TOL)
%

if nargin == 2
    tol = 0;
end

[~,~,info,~,~] = gklf(A,E,tol);
%

clc
disp(info);

disp('>> Right Kronecker structure');
lr = length(info.mr);
for i=1:lr
    disp([' ',num2str(info.nr(i)-info.mr(i)),' L[',num2str(i-1),',',num2str(i),'] blocks']);
end
disp(' ');
% 
disp('>> Infinite Jordan structure');
linf = length(info.minf);
for i=1:linf
    disp([' Jinf','[',num2str(info.minf(i)),']']);
end
disp(' ');
% 
disp('>> Finite Jordan structure');
if info.mf(1) ~= 0
    lf = length(info.mf);
    for i=1:lf
        disp([' J[',num2str(num2str(info.mf(i))),']']);
    end
end
disp(' ');
% 
disp('>> Left Kronecker structure');
ll = length(info.ml);
for i=ll:-1:1
    disp([' ',num2str(info.ml(i)-info.nl(i)),' LT[',num2str(ll-(i)),',',num2str(ll-(i-1)),'] blocks']);
end
disp(' ');

end
