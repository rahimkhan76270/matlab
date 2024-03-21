n=[7,12];
for i=1:2
    A=hilb(n(i));
    [Q1,R1]=cgsqr(A);
    [Q2,R2]=mgsqr(A);
    [Q3,R3]=qr(A,0);
    fprintf('n=%d\n',n(i));
    norm1=norm(Q1'*Q1-eye(n(i)),2);
    norm2=norm(Q2'*Q2-eye(n(i)),2);
    norm3=norm(A-Q3*R3,2);
    fprintf('||Q1^TQ1-I||=%e\n',norm1);
    fprintf('||Q2^TQ2-I||=%e\n',norm2);
    fprintf('||H-Q*R||=%e\n',norm3);
    fprintf('u*cond(A)=%e\n',eps*cond(A));
end
function [Q,R]=cgsqr(A)
    [~,n]=size(A);
    Q=A;
    R=zeros(n,n);
    for k=1:n
        R(1:k-1,k)=Q(:,1:k-1)'*A(:,k);
        Q(:,k)=A(:,k)-Q(:,1:k-1)*R(1:k-1,k);
        R(k,k)=norm(Q(:,k));
        Q(:,k)=Q(:,k)/R(k,k);
    end
end

function [Q,R]=mgsqr(A)
    [~,n]=size(A);
    Q=A;
    R=zeros(n,n);
    for k=1:n
        R(k,k)=norm(Q(:,k));
        Q(:,k)=Q(:,k)/R(k,k);
        R(k,k+1:n)=Q(:,k)'*Q(:,k+1:n);
        Q(:,k+1:n)=Q(:,k+1:n)-Q(:,k)*R(k,k+1:n);
    end
end