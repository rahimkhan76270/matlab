A=[4,-1,0;-2,6,1;0,1,7];
B=[3;9;6];
Ab=[A,B];
n=length(B);

for v=1:n-1
    Ab(v+1,:)=Ab(v+1,:)-(Ab(v+1,v)/Ab(v,v))*Ab(v,:);
end

x=zeros(n,1);
x(n)=Ab(n,end)/Ab(n,n);
for q=n-1:-1:1
    x(q)=(Ab(q,end)-Ab(q,q+1)*x(q+1))/Ab(q,q);
end

disp(x)