
n=8;
A=rand(n);
disp("Matrix");
disp(A);
disp("column max");
disp(max(A));
disp("row max");
disp(max(A,[],2));
disp("overall");
disp(max(max(A)));

disp("greater than 0.25");

for i=1:n
    for j=1:n
        if A(i,j)>0.25
            fprintf("%d %d \n",i,j);
        end
    end
end

