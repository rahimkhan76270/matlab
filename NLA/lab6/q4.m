year=1:11;
year=year';
price=[86.4;89.8;92.8;96.0;99.6;103.1;106.3;109.5;113.3;120.0;129.5];
[coeffs,residual] = solve_lsp(year,price);
disp(coeffs);
disp(residual);
fprintf('equation of least square line is y=%f+%fx\n',coeffs(1),coeffs(2));
fprintf('median house price in 2005 is %f\n',coeffs(1)+coeffs(2)*17);
fprintf('median house price in 2010 is %f\n',coeffs(1)+coeffs(2)*22);
function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end