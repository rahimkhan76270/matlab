A = [2 -1 0 0 -1 0 0 0;
    -1 3 -1 0 0 -1 0 0;
    0 -1 3 -1 0 0 -1 0;
    0 0 -1 2 0 0 0 -1;
    -1 0 0 0 2 -1 0 0;
    0 -1 0 0 -1 3 -1 0;
    0 0 -1 0 0 -1 3 -1;
    0 0 0 -1 0 0 -1 2];

B = [0 0 0 100 0 0 0 100];
n = length(B);

x = zeros(1, n, 'double');

for iter = 1:100
    for i = 1:n
        x(i) = B(i);
        for j = 1:n
            if i ~= j
                x(i) = x(i) - A(i, j) * x(j);
            end
        end
        x(i) = x(i) / A(i, i);
    end
end

disp(x);
