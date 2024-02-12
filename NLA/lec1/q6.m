n=16;
A=[-2,1,zeros(1,n-3),1];
B=toeplitz(A);
C=[1,2,3,4,5,6,7,8];
toeplitz(C)

format rational;

D=[1 1/2 1/3 1/4 1/5 1/6 1/7 1/8];
disp(toeplitz(D));