%% equation Ut=KUxx  k=1
% r=1 Dx=0.2 a=0 b=1
% T1=0 T2=1

format long e;
% initialization
r=1/6;
Nx=5;
Dx=0.2;
Dt=r*Dx*Dx;
Nt=round(1/Dt);
x=0:Dx:1;
u(:,1)=0.8*sin(pi*x);
T=0:Dt:1;

u(1,:)=0;
u(Nx+1,:)=0;
tic
for k=1:Nt
    for i=2:Nx
        u(i,k+1)=u(i,k)+r*(u(i-1,k)-2*u(i,k)+u(i+1,k));
    end
end
toc

disp(toc-tic)