A=[10^-16 1;1 1];
[L,U]=GENP(A);
E=L*U-A;
[L,U,p]=GEPP(A);
F=L*U-A(p,:);
disp(E)
disp(F)
disp(norm(E));
disp(norm(F));
disp(norm(E)/norm(A));
disp(norm(F)/norm(A));
%% solving Ax=b by GENP
A=[10^-16 1;1 1];
b=A*ones(2,1);
[L,U]=GENP(A);
y=Lower_triangular_solve(L,b);
xn=Upper_triangular_solve(U,y);
disp(xn);
err=norm(xn-ones(2,1))/norm(ones(2,1));
disp(err);
%% solving Ax=b by GEPP
A=[10^-16 1;1 1];
b=A*ones(2,1);
[L,U,p]=GEPP(A);
y=Lower_triangular_solve(L,b);
xp=Upper_triangular_solve(U,y);
disp(xp);
err=norm(xp-ones(2,1))/norm(ones(2,1));
disp(err);
%%
n=8;
H=hilb(n);
[L,U,P,d]=GEPP2(H);
disp(d);
disp(det(H))
%%
n=10;
H=hilb(n);
j=randn(n,1);
b=H*j;
x1=H\b;
r=H*x1-b;
disp([norm(r) norm(j-x1)])
n=[8 10 12];
for j=n
    H=hilb(j);
    [L,U]=GENP(H);
    E=L*U-H;
    err=norm(E)/norm(H);
    disp(err);
    x=randn(j,1);
    b=H*x;
    xx=U\(L\b);
    fprintf('relative error %d \n',norm(x-xx)/norm(x));
end
%%
function [L,U,P,d]=GEPP2(A)
    d=det(A);
    [L,U,P]=GEPP(A);
end
function [L, U] = GENP(A)
    [~, n] = size(A);
    for k = 1:n-1
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        j = k+1:n;
        A(j,j) = A(j,j)-A(j,k)*A(k,j);
    end
    L = eye(n,n)+ tril(A,-1);
    U = triu(A); 
end


function [L, U, p] = GEPP(A)
    [~, n] = size(A);
    p = (1:n)';
    for k = 1:n-1
        [~, m] = max( abs( A(k:n,k) ) );
        m = m+k-1;
        if (m ~=k)
            A([k m],:) = A([m k],:);
            p([k m]) = p([m k]);
        end
        if (A(k,k) ~= 0)
            A(k+1:n,k) = A(k+1:n,k)/A(k,k);
            j = k+1:n;
            A(j,j) = A(j,j)-A(j,k)*A(k,j);
        end
    end
    L = eye(n,n)+ tril(A,-1);
    U = triu(A); 
end
function vals=Upper_triangular_solve(A,b)
    N=length(A);
    vals=zeros(N,1);
    vals(N)=b(N)/A(N,N);
    for i=N-1:-1:1
        subt=0;
        for j=i+1:N
            subt=subt+A(i,j)*vals(j);
        end
        vals(i)=(b(i)-subt)/A(i,i);
    end
end

function vals=Lower_triangular_solve(A,b)
    N=length(A);
    vals=zeros(length(A),1);
    vals(1)=b(1)/A(1,1);
    for i=2:N
        subt=0;
        for j=1:i-1
            subt=subt+A(i,j)*vals(j);
        end
        vals(i)=(b(i)-subt)/A(i,i);
    end
end
