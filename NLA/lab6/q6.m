amount=[86;99.8;115.8;125;132.6;143.1;156.3;169.5];
year=1:8;
year=year';
[coeffs1,residual] = solve_lsp(year,amount);
fprintf('The coefficients are y= %f+%fx\n',coeffs1(1),coeffs1(2));
fprintf('amount of waste in 2005 is %f\n',coeffs1(1)+coeffs1(2)*9);
fprintf('amount of waste in 2010 is %f\n',coeffs1(1)+coeffs1(2)*10);
log_amount=log(amount);
[coeffs2,residual] = solve_lsp(year,log_amount);
fprintf('The coefficients are logy= %f+%fx\n',coeffs2(1),coeffs2(2));
fprintf('amount of waste in 2005 is %f\n',exp(coeffs2(1)+coeffs2(2)*9));
fprintf('amount of waste in 2010 is %f\n',exp(coeffs2(1)+coeffs2(2)*10));
plot(year,amount,'o',year,exp(coeffs2(1)+coeffs2(2)*year));
hold on;
plot(year,amount,'o',year,coeffs1(1)+coeffs1(2)*year);
hold off;
function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end