a=0;b=1;
C=0.85;
alpha=0.5;
Dx=0.1;
Dt=Dx*C/alpha;
Nx=round((b-a)/Dx);
Nt=round(1/Dt);
u=zeros(Nt+1,Nx+1);
for i=1:Nt+1
    u(i,Nx+1)=2-(i-1)*Dt;
end
x=a:Dx:b;
u(1,:)=2*x;
for i=2:Nt
    for j=2:Nx
        if mod(i,2)==0
            u(i,j)=u(i-1,j)-C*(u(i-1,j+1)-u(i-1,j));
        else
             u(i,j)=u(i-1,j)-0.5*(0.5*(u(i-1,j+1)-u(i-1,j-1))+0.5*C*(u(i-1,j+1)-u(i-1,j-1)));
        end 
    end
end
vals=exact(a,b,1,Dx,Dt);
disp(norm(u-vals,2));
function vals=exact(a,b,t,Dx,Dt)
    Nx=round((b-a)/Dx);
    Nt=round(t/Dt);
    vals=zeros(Nt+1,Nx+1);
    for i=1:Nt+1
        for j=1:Nx+1
            vals(i,j)=2*(j-1)*Dx-(i-1)*Dt;
        end
    end
end