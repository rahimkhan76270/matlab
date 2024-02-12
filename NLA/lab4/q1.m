

for n=[8,10,12]
    H=hilb(n);
    HI=invhilb(n);
    x=ones(n,1);
    b=H*x;
    x1=H\b;
    x2=HI*b;
    disp('x1')
    disp(x1);
    disp('x2');
    disp(x2);
    
    disp('for method 1')
    eta=norm(H*x1-x)/(norm(H)*norm(x));
    err=norm(x1-x)/norm(x);
    cond_num1=cond(H);
    disp('eta');
    disp(eta);
    disp('conditon number')
    disp(cond_num1);
    disp('error')
    disp(err);
    
    eta=norm(H*x2-x)/(norm(H)*norm(x));
    err=norm(x2-x)/norm(x);
    cond_num1=cond(H);
    disp('eta');
    disp(eta);
    disp('conditon number')
    disp(cond_num1);
    disp('error')
    disp(err);
end


n=10;
H=hilb(n);
x=rand(n,1);
b=H*x;
x1=H\b;
r=H*x1-b;
disp([norm(r,inf),norm(x-x1,inf)]);


n=10;
A=diag(ones(n,1)*2)+diag(ones(n-1,1)*(-1),-1)+diag(ones(n-1,1)*(-1),1);
A(1,1)=1;
[L,U]=GENP(A);
disp(L);
disp(U);
D=diag(diag(A));
A1=L*D*L';
disp(A1);
ch=chol(A);
disp(ch);



% Matlab program that implements growth of cond(H)

% generate Hilbert matrices and compute cond number with 2-norm

 

N = 50; % maximum size of a matrix

condofH = []; % conditional number of Hilbert Matrix

N_it = zeros(1, N);

 

% compute the cond number of Hn

for n = 1:N

    Hn = hilb(n);

    N_it(n) = n;

    condofH = [condofH cond(Hn, 2)];

end

 

x = 1:50;

y = (1 + sqrt(2)).^(4 * x) ./ (sqrt(x));

 

% plot

plot(N_it, log(condofH), x, log(y));

hold on; % to overlay plots

plot(N_it, log(condofH), 'b', x, log(y), 'r');

hold off; % release the hold

 

% plot labels

title('Conditional Number growth of Hilbert Matrix: Theoretical vs Matlab');

xlabel('N', 'fontsize', 16);

ylabel('log(cond(Hn))', 'fontsize', 16);

lgd = legend({'MatLab', 'Theoretical'}, 'Location', 'northwest');

legend('show');



function [L, U] = GENP(A)
    [~, n] = size(A);
    for k = 1:n-1
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        j = k+1:n;
        A(j,j) = A(j,j)-A(j,k)*A(k,j);
    end
    L = eye(n,n)+ tril(A,-1);
    U = triu(A); 
end
