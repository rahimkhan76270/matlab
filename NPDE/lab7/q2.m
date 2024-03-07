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
    a=ones(Nx-1,1)*r;
    b=ones(Nx-1,1)*(1+2*r);
    c=ones(Nx-1,1)*r;
    a(1)=0;
    c(Nx-1)=0;
    u(j,2:Nx)=thomas_algorithm(a,b,c,d);
end
exact_solution=exact_value(Dx,Dt);
err=norm(exact_solution-u,2);
disp('error');
disp(err);
surf(u);


function val=exact_value(Dx,Dt)
    a=0;b=4;
    Nx=round((b-a)/Dx);
    Nt=round(1/Dt);
    val=zeros(Nt+1,Nx+1);
    exact_sol=@(x,t)(exp(-t)*sin(pi*x/2))+exp(-t/4)*sin(pi*x/4);
    for i=1:Nt+1
       for j=1:Nx+1 
           val(i,j)=exact_sol((j-1)*Dx,(i-1)*Dt);
       end
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
