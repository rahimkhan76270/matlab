%% q1
h=[0.282; 0.564; 0.752; 0.940];
d=[0.752; 1.102; 1.248; 1.410];
[coeffs1,residual1] = solve_lsp(h,d);
disp(coeffs1);
hmax=max(h);
hmin=min(h);
h1=linspace(hmin,hmax,100);
plot(h,d,'b');
hold on;
plot(h1,coeffs1(1)+coeffs1(2)*h1,'r');
hold off;

function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end