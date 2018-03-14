clear all
close all

im1 = imread('Slc1-06.gif');
im2 = imread('Slc1-07.gif');

Ixfilt = [-1 2 -1; 0 0 0; -1 2 1];
Ix = conv2(double(im1), double(Ixfilt));
Iy = conv2(double(im1), double(Ixfilt'));
It = double(im2-im1);
a = .1;
u_o = zeros(74,74);
v_o = zeros(74,74);
H = [1/12 1/6 1/12; 1/6 -1 1/6; 1/12 1/6 1/12];

diff = 1;
tol = .01;
k = 1;

% jacobi method
while and(diff > tol, k < 100)
    u_avg = conv2(double(u_o), H);
    v_avg = conv2(double(v_o), H);
    u = nan(74,74);
    v = nan(74,74);
    for i = 1:74
        for j = 1:74
            u(i,j) = u_avg(i,j) - (Ix(i,j)*(Ix(i,j)*u_avg(i,j) + Iy(i,j)*v_avg(i,j) + It(i,j)))...
                /(a^2 + Ix(i,j)^2 + Iy(i,j)^2);
            v(i,j) = v_avg(i,j) - (Iy(i,j)*(Ix(i,j)*u_avg(i,j) + Iy(i,j)*v_avg(i,j) + It(i,j)))...
                /(a^2 + Ix(i,j)^2 + Iy(i,j)^2);
        end
    end
    diff_v(k) = 0;
    for m = 1:74
        for n = 1:74
            vec1 = [u(m,n) v(m,n)];
            vec2 = [u_o(m,n) v_o(m,n)];
            diff_v(k) = diff_v(k) + acos(min(1,max(-1,vec1(:).'*vec2(:)/norm(vec1)/norm(vec2))));
        end
    end
    diff_v(k) = diff_v(k)/(74*74);
    diff = diff_v(k);
    v_o = v;
    u_o = u;
    k = k+1;
end

[x,y] = meshgrid(1:74,1:74);

figure(1)
hold on;
quiver(x,y,u,v);
xlim([1 74]);
ylim([1 74]);
title('Velocity Field Solved with Jacobi Method');

figure(2)
hold on;
plot(diff_v);
title('Average Per-Pixel Difference in Direction of Velocity Field (Jacobi)')
xlabel('Iteration')
ylabel('Average Difference')

u_j = u;
v_j = v;

% gauss-sidel method
u = zeros(74,74);
v = zeros(74,74);
diff=1;
diff_v = 0;
k = 1;

while and(diff > tol, k < 100)
    u_o = u;
    v_o = v;
    for i = 1:74
        for j = 1:74
            u_avg = conv2(double(u), H);
            v_avg = conv2(double(v), H);
   	        u(i,j) = u_avg(i,j) - (Ix(i,j)*(Ix(i,j)*u_avg(i,j) + Iy(i,j)*v_avg(i,j) + It(i,j)))...
                /(a^2 + Ix(i,j)^2 + Iy(i,j)^2);
            v(i,j) = v_avg(i,j) - (Iy(i,j)*(Ix(i,j)*u_avg(i,j) + Iy(i,j)*v_avg(i,j) + It(i,j)))...
                /(a^2 + Ix(i,j)^2 + Iy(i,j)^2);
        end
    end
    diff_v(k) = 0;
    for m = 1:74
        for n = 1:74
            vec1 = [u(m,n) v(m,n)];
            vec2 = [u_o(m,n) v_o(m,n)];
            diff_v(k) = diff_v(k) + acos(min(1,max(-1,vec1(:).'*vec2(:)/norm(vec1)/norm(vec2))));
        end
    end
    diff_v(k) = diff_v(k)/(74*74);
    diff = diff_v(k);
    k = k+1;
end

figure(3)
hold on;
quiver(x,y,u,v);
xlim([1 74]);
ylim([1 74]);
title('Velocity Field Solved with Gauss-Sidel Method');

figure(4)
hold on;
plot(diff_v);
title('Average Per-Pixel Difference in Direction of Velocity Field (Gauss-Sidel)')
xlabel('Iteration')
ylabel('Average Difference')
