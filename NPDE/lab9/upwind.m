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
        %1st order upwind unstable
        % u(i+1,j)=u(i,j)-0.5*C*(3*u(i,j)-4*u(i,j-1)-u(i,j-2)); 
        %2nd order upwind stable
        u(i+1,j)=u(i,j)-C*(u(i,j+1)-u(i,j-1))-0.5*C*(1-C)*(u(i,j)-2*u(i,j-1)+u(i,j-2)); 
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