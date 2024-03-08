temp=[72;79;88;96];
time=[0;20;40;60];
[coeffs,residual] = solve_lsp(time,temp);
disp((165-coeffs(1))/coeffs(2));
function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end