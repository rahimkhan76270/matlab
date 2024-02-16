k=1;Dx=0.1;r=1/2;
Nx=round(1/Dx);
Dt=r*Dx*Dx/k;
Nt=round(1/Dt);
x=0:Dx:1;
u=zeros(Nt+1,Nx+1);
u2=zeros(Nt+1,Nx+1);
u(1,:)=x.*(1-x);
for j=2:Nt+1
    d=zeros(Nx+1,1);
    d(1)=r*u(j-1,2)+(1-r)*u(j,1)-2*r*Dx+Dt*cos(0);
    d(Nx+1)=u(j-1,Nx+1)+(1-r)*u(j,Nx)+2*r*Dx+Dt*cos((Nx+1)*Dx);
    for i=2:Nx
        d(i-1)=r*u(j-1,i+1)+(1-2*r)*u(j-1,i)+r*u(j-1,i-1)+Dt*cos((i-1)*Dx)-2*r*Dx;
    end
    A=diag(ones(Nx+1,1)*(1+2*r))+diag(ones(Nx,1)*(-r),-1)+diag(ones(Nx,1)*(-r),1);
    A(1,1)=1+r;
    A(1,2)=-r;
    A(Nx+1,Nx+1)=1+r;
    A(Nx+1,Nx)=-r;
    u(j,1:Nx+1)=tridiagonal_solve(A,d);
    b=ones(Nx+1,1)*(1+2*r);
    a=ones(Nx+1,1)*(r);
    c=ones(Nx+1,1)*(r);
    b(1)=1+r;
    b(Nx+1)=1+r;
    c(Nx+1)=0;
    a(1)=0;
    u2(j,1:Nx+1)=thomas_algorithm(a,b,c,d);
end

disp(norm(u-u2,2));
function vals=tridiagonal_solve(A,b)
    Aug=[A,b];
    N=length(A);
    for i=2:N
        Aug(i,:)=Aug(i,:)-Aug(i-1,:)*Aug(i,i-1)/Aug(i-1,i-1);
    end
    vals=Upper_triangular_solve(Aug(1:N,1:N),Aug(:,N+1));
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