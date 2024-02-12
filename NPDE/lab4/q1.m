k=4/(pi*pi);Dx=0.1;r=1/2;
Nx=round(4/Dx);
Dt=r*Dx*Dx/k;
Nt=round(1/Dt);
x=0:Dx:4;
u=zeros(Nt+1,Nx+1);
u(1,:)=function1(x);
u(:,1)=0;
u(:,Nx+1)=0;
for j=2:Nt+1
    b=zeros(Nx-1,1);
    for i=2:Nx
        b(i-1)=r*u(j-1,i+1)+(1-2*r)*u(j-1,i)+r*u(j-1,i-1);
    end
    A=diag(ones(Nx-1,1)*(1+2*r))+diag(ones(Nx-2,1)*(-r),-1)+diag(ones(Nx-2,1)*(-r),1);
    u(j,2:Nx)=tridiagonal_solve(A,b);
end
exact_solution=function2(Dx,Dt);
err=norm(exact_solution-u,2);
disp('error');
disp(err);
surf(u);
% surf(exact_solution);
function vals=tridiagonal_solve(A,b)
    Aug=[A,b];
    N=length(A);
    for i=2:N
        Aug(i,:)=Aug(i,:)-Aug(i-1,:)*Aug(i,i-1)/Aug(i-1,i-1);
    end
    vals=Upper_triangular_solve(Aug(1:N,1:N),Aug(:,N+1));
end

function val=function1(x)
    val=sin(pi*x/4).*(1+2*cos(pi*x/4));
end
function val=function2(Dx,Dt)
    m=round(4/Dx);
    n=round(1/Dt);
    val=zeros(n+1,m+1);
    for x=1:m+1
        for t=1:n+1
            val(t,x)=(exp(-t*Dt).*sin(pi*(x-1)*Dx/2)+exp(-t*Dt/4).*sin(pi*(x-1)*Dx/4));
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
