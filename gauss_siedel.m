A=[4,-1,-1;-2,6,1;-1,1,7];
B=[3,9,-6];

values=zeros(10,3,'double');

prev=[0,0,0];
next=[1,1,1];
for i=1:10
    x=(B(1)-A(1,2)*y-A(1,3)*z)/A(1,1);
    y=(B(2)-A(2,1)*x-A(2,3)*z)/A(2,2);
    z=(B(3)-A(3,1)*x-A(3,2)*y)/A(3,3);
    values(i,1)=x;
    values(i,2)=y;
    values(i,3)=z;
end
% disp(values);

for i=1:10
     for j=1:3
         prev(j)=B(j);
         for k=1:3
             if j~=k
                 prev(j)=prev(j)-A(j,k)*prev(k);
             end
         end
         prev(j)=prev(j)/A(j,j);
     end
end
disp(prev);