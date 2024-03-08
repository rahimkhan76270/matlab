year=0:5;
year=year';
population=[76;106;132;181;227;282];
log_population=log(population);
[coeffs,residual]=solve_lsp(year,log_population);
fprintf('the population in 2020 is %f\n',exp(coeffs(1)+coeffs(2)*6));
fprintf('the population in 2050 is %f\n',exp(coeffs(1)+coeffs(2)*7.5));
function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end