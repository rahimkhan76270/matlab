x=[12;14;17;21;26;30];
y=[60;70;90;100;100;120];
[coeffs,residual] = solve_lsp(x,y);
fprintf('The coefficients are %f and %f\n',coeffs(1),coeffs(2));
fprintf('The residual is %f\n',residual);
fprintf('The line is y=%f+%fx\n',coeffs(1),coeffs(2));
fprintf('profit at 50000 = %f\n',coeffs(1)+coeffs(2)*50);
fprintf('profit at 100000 = %f\n',coeffs(1)+coeffs(2)*100);
xmax=max(x);
xmin=min(x);
xplot=linspace(xmin,xmax,100);
plot(x,y,'+');
hold on;
plot(xplot,coeffs(1)+coeffs(2)*xplot);


function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end