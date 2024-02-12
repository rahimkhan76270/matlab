
A=[3 ,0, 0, 0, 0, -1, -1, -1;
0 2 0 -1 0 0 -1 0;
0 0 3 -1 0 0 -1 -1;
0 -1 -1 2 0 0 0 0;
0 0 0 0 2 -1 0 -1;
-1 0 0 0 -1 2 0 0;
-1 -1 -1 0 0 0 3 0;
-1 0 -1 0 -1 0 0 3];

% A=diag([1,1,2,])

B=[0 100 0 100 0 0 0 0];
n=length(B);

prev=zeros(1,n,"double");
for i=1:10
     for j=1:n
         prev(j)=B(j);
         for k=1:n
             if j~=k
                 prev(j)=prev(j)-A(j,k)*prev(k);
             end
         end
         prev(j)=prev(j)/A(j,j);
     end
end

disp(prev);
%%
% d=det(A)
% C=inv(A);
% Y=inv(A).*B;
%%