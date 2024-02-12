k=1;Dx=0.1;r=1/2;
Nx=round(1/Dx);
Dt=r*Dx*Dx/k;
Nt=round(1/Dt);
x=0:Dx:1;
u=zeros(Nt+1,Nx+1);
u(1,:)=function1(x);
u(:,1)=0;
u(:,Nx+1)=0;
for j=2:Nt+1
    d=zeros(Nx-1,1);
    for i=2:Nx
        d(i-1)=r*u(j-1,i+1)+(1-2*r)*u(j-1,i)+r*u(j-1,i-1);
    end
    % A=diag(ones(Nx-1,1)*(1+2*r))+diag(ones(Nx-2,1)*(-r),-1)+diag(ones(Nx-2,1)*(-r),1);
    a=ones(Nx-1,1)*r;
    b=ones(Nx-1,1)*(1+2*r);
    c=ones(Nx-1,1)*r;
    a(1)=0;
    c(Nx-1)=0;
    % u(j,2:Nx)=tridiagonal_solve(A,d);
    u(j,2:Nx)=thomas_algorithm(a,b,c,d);
end
exact_solution=function2(Dx,Dt);
err=norm(exact_solution-u,2);
disp('error');
disp(err);
function vals=tridiagonal_solve(A,b)
    Aug=[A,b];
    N=length(A);
    for i=2:N
        Aug(i,:)=Aug(i,:)-Aug(i-1,:)*Aug(i,i-1)/Aug(i-1,i-1);
    end
    vals=Upper_triangular_solve(Aug(1:N,1:N),Aug(:,N+1));
end

function val=function1(x)
    val=sin(pi*x).*(1+2*cos(pi*x));
end
function val=function2(Dx,Dt)
    m=round(1/Dx);
    n=round(1/Dt);
    val=zeros(n+1,m+1);
    for x=1:m+1
        for t=1:n+1
            val(t,x)=(exp(-4*pi*pi*(t-1)*Dt).*sin(2*pi*(x-1)*Dx)+exp(-pi*pi*(t-1)*Dt).*sin(pi*(x-1)*Dx));
        end
    end
end

function vals=Upper_triangular_solve(A,b)
    N=length(A);
    vals=zeros(N,1);
    vals(N)=b(N)/A(N,N);
    for i=N-1:-1:1
        vals(i)=(b(i)-A(i,i+1)*vals(i+1))/A(i,i);
    end
end

function val=thomas_algorithm(a,b,c,d)
% [b -c 0]
% [-a b -c]
% [0 -a b]
% all the negative terms are considered in the calculation just give the +ve values of a and c 
    N=length(a);
    val=zeros(N,1);
    alpha=zeros(N,1);
    s=zeros(N,1);
    alpha(1)=b(1);
    s(1)=d(1);
    for j=2:N
       alpha(j)=b(j)-a(j)*c(j-1)/alpha(j-1);
    end
    for j=2:N
       s(j)=d(j)+a(j)*s(j-1)/alpha(j-1);
    end
    val(N)=s(N)/alpha(N);
    for j=N-1:-1:1
       val(j)=(s(j)+c(j)*val(j+1))/alpha(j); 
    end
end
