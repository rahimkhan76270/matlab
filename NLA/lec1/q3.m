n=5;
A=magic(n);
disp(A);
p=randperm(n);
disp(p);

q=randperm(n);
disp(q);
A=A(p,q);
disp(A);


disp(sum(A))
disp(sum(transpose(A)))
disp(transpose(sum(transpose(A))))
disp(sum(diag(A)))
disp(sum(diag(flipud(A))))

disp(rank(A))



A=magic(4);
disp(A);

disp(null(A));
disp(null(A,'rational'));

disp(null(sym(A)))

rref(A)
rank(A)