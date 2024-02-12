M=magic(5);
disp(M);
disp("column sums")
disp(sum(M));
disp("row sum");
disp(sum(M,2));
disp("diagonal sum")
disp(sum(diag(M)));
X=flip(M);
disp("anti diagonal sum");
disp(sum(diag(X)));
