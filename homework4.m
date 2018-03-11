x = linspace(-10, 10)
o = .7
c = 1
w = 3

g = exp((-1/(2*o^2))*(x-c).^2)
g_w = exp((-1/(2*o^2))*(x-c).^2).*w

X = [1;-1]
y = [1; 1]
N = size(y, 1)
c_1 = 1
c_2 = linspace(-10, 10)
J = zeros(size(c_2))
w_1 = 1
w_2 = 1

for i = 1:size(c, 2)
    J(i) = (1/N)*sum((y-(exp((-1/(2*o^2))*(X-c_1).^2).*w_1 +exp((-1/(2*o^2))*(X-c_2(i)).^2).*w_2)).^2)
end

figure(1)
hold on
plot(x, g)
plot(x, g_w)
xlabel('x')
ylabel('g(x)')
legend('w = 1', 'w = 3')
title('g(x).w for sigma=1, c=1')

figure(2)
hold on
plot(c_2, J)
xlabel('c2 (center for second basis function)')
ylabel('J')
title('Loss function J for sigma=0.7, w1=1, w2=1, c1=1')