[U,X]=qr(randn(80));
[V,X]=qr(randn(80));
S=diag(2^(-1:-1:-80));
A=U*S*V';
[QC,RC]=cgsqr(A);
[QM,RM]=mgsqr(A);
[Q,R]=qr(A);
fprintf('||Q-QC||=%e\n',norm(QC'*QC-eye(80)));
fprintf('||Q-QM||=%e\n',norm(QM'*QM-eye(80)));
fprintf('||Q-Q||=%e\n',norm(Q'*Q-eye(80)));
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