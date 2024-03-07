h=15;
w=10;
T_top=100;
T_left=0;
T_right=0;
T_bottom=0;
Dx=0.5;
Dy=0.5;
beta=Dx/Dy;
Nx=round(w/Dx);
Ny=round(h/Dy);
u=zeros(Ny+1,Nx+1,'double');
u(Ny+1,:)=T_top;
u(1,:)=T_bottom;
u(:,1)=T_left;
u(:,Nx+1)=T_right;
sz=Nx*Ny;
A=diag(ones(sz,1)*(-2-2*beta*beta))+diag(ones(Nx*Ny-1,1)*(beta*beta),1)+diag(ones(Nx*Ny-1,1)*beta,-1);%+diag(ones(Nx*Ny-2*Ny-4),Nx+1)+diag(ones(Nx*Ny-2*Ny-4),-Nx-1);

for i=1:Ny-1
   A(i*Nx,i*Nx+1)=0;
   A(i*Nx+1,i*Nx)=0;
end
b=zeros(Nx*Ny,1);
for i=1:Ny
   b(i*Nx)=-beta*beta*T_top; 
end


function vals=actual_value(h,w,Dx,Dy)
    Nx=round(w/Dx);
    Ny=round(h/Dy);
    vals=zeros(Ny+1,Nx+1);
    for i=1:Ny+1
       for j=1:Nx+1
            n=20;
            s=0;
            for k=1:n
                s=s+(sinh(k*pi*(j-1)*Dy/w))*(sin(k*pi*(j-1)*Dx/w))/(k*sinh(k*pi*h/w));
            end
            vals(i,j)=400*s/pi;
       end
    end
end
