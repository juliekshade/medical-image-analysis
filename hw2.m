clear all
close all

im = imread('xid-23155723_2.bmp');
im = im(:, :, 1);

a = 20
b = 200

L = double(max(max(im)));

figure(1);
h_orig = histogram(im, L);

h_orig_counts = h_orig.Values;
total_px = size(im, 1)*size(im, 2)
h_x = h_orig_counts./total_px;

for i=1:L
    H_x(i) = sum(h_x(1:i));
    h_z(i) = b + a*i;
end

h_z = h_z./(sum(h_z));

for i=1:L
    H_z(i) = sum(h_z(1:i));
end

im_new = im;
lookuptable = zeros(L, 1);
final_hist = zeros(L, 1);

for i=1:L
    diff = abs(H_x(i) - H_z);
    [M, I] = min(diff);
    lookuptable(i) = I;
    final_hist(lookuptable(i)) = final_hist(lookuptable(i)) + h_orig_counts(i);
    im_new(im == i) = lookuptable(i);
end

figure(2)
hold on
histogram(im_new, L);
plot([0:1:L-1], h_z*total_px);

figure(3)
imshow(im_new);