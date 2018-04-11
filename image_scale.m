close all
clear all

% read both images
im1 = imread('img1.gif');
im2 = imread('img2.gif');

% normalize images and set inital values of scale and MSE
im1 = im1./max(max(im1));
im2 = im2./max(max(im2));
s = 1;
s_old = 0;
mse = 0;
i = 1;

while abs(s-s_old) > .001
im1_sc = imresize(im1,s); % scale image
d = size(im1_sc,1); 
if mod(d,2)==0 % odd number -> make D even number of pixels
    F = griddedInterpolant(double(im1_sc));
    [sx,sy,sz] = size(im1_sc);
    [xq,yq] = ndgrid(0:(d-1)/d:sx, 0:(d-1)/d:sy);
    im1_sc = uint8(F(xq,yq)); % interpolate image so can find px around center
    d=d+1;
end
% compare overlapping regions
if d >= 256
    im1_sc_d = im1_sc(d/2-127:d/2+128, d/2-127:d/2+128);
    msenew = sum(sum((im2-im1_sc_d).^2))/(size(im2,1)^2);
else
    im2_sc = im2(129-d/2:256-(128-d/2), 129-d/2:256-(128-d/2));
    im1_sc = im1_sc(129-d/2:256-(128-d/2), 129-d/2:256-(128-d/2));
    msenew = sum(sum((im2_sc-im1_sc).^2))/(size(im2_sc,1)^2);
end
s_old = s
s = s_old - (2/i)*(-msenew+mse); % update s (by less each iteration)
mse = msenew;
mse_all(i) = mse; 
i =i+1;
end

figure(1)
plot(mse_all)
title('Convergence of gradient descent');
xlabel('Iteration Number')
ylabel('Mean Squared Error')