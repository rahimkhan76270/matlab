accounts = 100 + (100000-100) * rand(50000,1);
accounts = floor(100 * accounts)/100;   
ill_acc=0;

d=0;
while ill_acc<1000000
    accounts=accounts*(1+5/(365*100));
    truncated=floor(100*accounts)/100;
    pennies=accounts-truncated;
    accounts=truncated;
    ill_acc=ill_acc+sum(pennies);
    ill_acc=ill_acc*(1+5/(365*100));
    d=d+1;
end
disp(d)
