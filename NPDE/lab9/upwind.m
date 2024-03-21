% using left boundary u(0,t)=0;
l=1;
C=0.95;
alpha=0.1;
Dx=0.1;
Dt=Dx*C/alpha;
Nx=round(l/Dx);
Nt=round(1/Dt);
u=zeros(Nt+1,Nx+1);
for j=1:Nx+1
    if (j-1)*Dx>3/8 && (j-1)*Dx <5/8
        u(1,j)=1;
    end
end

for i=2:Nt+1
    for j=2:Nx+1
        u(i,j)=u(i-1,j)-C*(u(i-1,j)-u(i-1,j-1));   
    end
end
vals=exact(0,l,1,Dx,Dt,alpha);
disp(norm(u-vals,2));
function vals=exact(a,b,t,Dx,Dt,alpha)
    Nx=round((b-a)/Dx);
    Nt=round(t/Dt);
    vals=zeros(Nt+1,Nx+1);
    for i=1:Nt+1
        for j=2:Nx+1
            if (j-1)*Dx-alpha*(i-1)*Dt >3/8 && (j-1)*Dx-alpha*(i-1)*Dt <5/8
                vals(i,j)=1;
            end
        end
    end
end