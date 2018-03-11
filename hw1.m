clear all
close all


im = analyze75read('KKI2009-01-MPRAGE_stripped.img');
D = double(im);
V = reshape(D, [], 1);
h = histcounts(V, 1:50000:max(V))
i = 1:50000:max(V)
hold on
histogram(V,1:50000:max(V))
xlabel('Pixel Value')
ylabel('Frequency')


GM_mask = im;
GM_mask(im < i(6)) = 0;
GM_mask(im > i(14)) = 0;
GM_mask = (logical(GM_mask));
imshow(GM_mask(:, :, 155));
M2 = reshape(GM_mask(:, 60, :), [256 256]);
T0 = maketform('affine',[0 -1; 1 0; 0 0]);
R2 = makeresampler({'cubic','nearest'},'fill');
M3 = imtransform(M2,T0,R2);  
figure, imshow(M3);


WM_mask = im;
WM_mask(im <= i(14)) = 0;
WM_mask = (logical(WM_mask));
imshow(WM_mask(:, :, 155));
M2 = reshape(WM_mask(:, 60, :), [256 256]);
T0 = maketform('affine',[0 -1; 1 0; 0 0]);
R2 = makeresampler({'cubic','nearest'},'fill');
M3 = imtransform(M2,T0,R2);  
figure, imshow(M3);

gm_vox = sum(reshape(GM_mask, [], 1))
wm_vox = sum(reshape(WM_mask, [], 1))