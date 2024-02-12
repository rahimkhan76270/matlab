a=ones(5,1);
a(1)=0;
a(5)=2;
b=ones(5,1)*4;
c=[1;1;1;1;0];
d=[0.4;0.8;1.2;1.6;1.6];
disp(thomas_algorithm(a,b,c,d));

function val=thomas_algorithm(a,b,c,d)
% [b -c 0]
% [-a b -c]
% [0 -a b]
% all the negative terms are considered in the calculation just give the +ve values of a and c 

    N=length(a);
    val=zeros(N,1);
    alpha=zeros(N,1);
    s=zeros(N,1);
    alpha(1)=b(1);
    s(1)=d(1);
    for j=2:N
       alpha(j)=b(j)-a(j)*c(j-1)/alpha(j-1);
    end
    for j=2:N
       s(j)=d(j)+a(j)*s(j-1)/alpha(j-1);
    end
    disp(alpha);
    disp(s);
    val(N)=s(N)/alpha(N);
    for j=N-1:-1:1
       val(j)=(s(j)+c(j)*val(j+1))/alpha(j); 
    end
end