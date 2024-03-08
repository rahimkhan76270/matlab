%% 5.5.1 a
t1=[-2;0;1;3];
y1=[0;1;2;5];
t2=[1;2;3;4;5];
y2=[1;0;-2;-3;-3];
t3=[-2;-1;0;1;2];
y3=[-5;-3;-2;0;3];
[coeffs1,residual1] = solve_lsp(t1,y1);
[coeffs2,residual2] = solve_lsp(t2,y2);
[coeffs3,residual3] = solve_lsp(t3,y3);
figure;
subplot(2,2,1);
plot(t1,y1,'b');
hold on;
plot(t1,coeffs1(1)+coeffs1(2)*t1,'r');
hold off;
title('Graph 1');

subplot(2,2,2);
plot(t2,y2,'b');
hold on;
plot(t2,coeffs2(1)+coeffs2(2)*t2,'r');
hold off;
title('Graph 2');

subplot(2,2,3);
plot(t3,y3,'b');
hold on;
plot(t3,coeffs3(1)+coeffs3(2)*t3,'r');
hold off;
title('Graph 3');
%%

function [coeffs,residual] = solve_lsp(x,y)
    A=[ones(size(x)),x];
    coeffs = A\y;
    residual = norm(A*coeffs-y,2);
end