fitnessfcn = @(x)[norm(x)^2,0.5*norm(x(:)-[2;-1])^2+2];

A = [1,1];
b = 1/2;

rng default % For reproducibility
x = gamultiobj(fitnessfcn,2,A,b);

figure;
plot(x(:,1),x(:,2),'ko')
t = linspace(-1/2,2);
y = 1/2 - t;
hold on
plot(t,y,'b--')
xlabel('x(1)')
ylabel('x(4)')
title('Pareto Points in Parameter Space')
hold off