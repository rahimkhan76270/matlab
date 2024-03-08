time=0:7;
time=time';
mg=[100;82.7;68.3;56.5;46.9;38.6;31.9;26.4];
log_mg=log(mg);
[coeffs,residual]=solve_lsp(time,log_mg);
fprintf('amount after 10 days %f\n',exp(coeffs(1)+coeffs(2)*10));
fprintf('time for safe disposal %f days\n',(log(0.01)-coeffs(1))/coeffs(2));
function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end