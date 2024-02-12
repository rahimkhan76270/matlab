%% (a)
f1=@(x) 1/(sqrt(1+x)+sqrt(x));
x=10000;
val=f1(x);
disp(round(val,5,"decimals"))
%% (b)
f2=@(c,d)(2*sin((c-d)/2)*cos((c+d)/2));
x=10^(-50);
y=10^(-51);
disp(f2(x,y));

%% (c)

f3=@(x)(tan(x/2));
x=1;
val=f3(x);
while x>0
    val=f3(x);
    x=x/10;
end
disp(val);