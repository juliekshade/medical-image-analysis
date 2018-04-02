clear all;
close all;

x = [0 0; 3 0; 0 2];
y = [0 -1; 2 0; 0 2];

x_bar = mean(x);
y_bar = mean(y);
X = bsxfun(@minus, x, x_bar)';
Y = bsxfun(@minus, y, y_bar)';

A = Y*X';
[U,S,V] = svd(A);
R = U*diag([1, det(U*V')])*V';
t = y_bar' - R*x_bar';

for i = 1:3
    y_pred(i,:) = [R*x(i,:)' + t]';
end

figure(1)
scatter(x(:,1), x(:,2), 40, 'filled')
hold on;
scatter(y(:,1), y(:,2), 40)
title('Landmarks before registration');
legend('x_i', 'y_i')

figure(2)
scatter(y_pred(:,1), y_pred(:,2), 40, 'filled')
hold on;
scatter(y(:,1), y(:,2), 40)
title('Landmarks after registration')
legend('R*x_i + t', 'y_i')

err = 0;
for i = 1:3
    err = err + norm(y(i,:)-x(i,:))^2;
end

posterr = 0;
for i = 1:3
    posterr = posterr + norm(y(i,:)-y_pred(i,:))^2;
end