clear all
close all

%load all images
indiv = 3;
imgs = nan(192,168,64);
for i = 1:64
    imgs(:,:,i) = loadimage(indiv,i);
end

a = reshape(imgs(:,:,:),[size(imgs,1)*size(imgs,2),64]);
[mu,Ud,Y,sigma] = pca(a,2);

eigenarray = [reshape(((mu-min(mu))./max(mu-min(mu))), [size(imgs,1),size(imgs,2)]),...
    reshape((Ud(:,1)-min(Ud(:,1)))./max(Ud(:,1)-min(Ud(:,1))),[size(imgs,1),size(imgs,2)]),...
    reshape((Ud(:,2)-min(Ud(:,2)))./max(Ud(:,2)-min(Ud(:,2))),[size(imgs,1),size(imgs,2)])];

% plot mean face u and first 2 eigenfaces (shifted to fix negative values)
figure(1)
montage(eigenarray,'size',[1 NaN]);
title('Mean Image, First Eigenimage, Second Eigenimage');

o1 = -sigma(1):.2*sigma(1):sigma(1);
o2 = -sigma(2):.2*sigma(2):sigma(2);
p1 = mu + bsxfun(@times, o1,Ud(:,1));
p2 = mu + bsxfun(@times, o2,Ud(:,2));

p1array = [];
p2array = [];
for i = 1:size(p1,2)
    p1array = [p1array,reshape((p1(:,i)-min(p1(:,i)))./max(p1(:,i)-min(p1(:,i))),...
        [size(imgs,1),size(imgs,2)])];
    p2array = [p2array,reshape((p2(:,i)-min(p2(:,i)))./max(p2(:,i)-min(p2(:,i))),...
        [size(imgs,1),size(imgs,2)])];
end

% plot mu + y*u for y1 and y2
figure(2)
montage(p1array,'size',[1 NaN]);
title('µ + y1u1 for y1 = o1 : 0.2o1 : o1')

figure(3)
montage(p2array,'size',[1 NaN]);
title('µ + y2u2 for y2 = o2 : 0.2o2 : o2')

function [mu,Ud,Y,sigma]=pca(X,d)
    mu = mean(X,2);
    X_shift = bsxfun(@minus, X, mu);
    [U,S,V] = svd(X_shift,0);
    Ud = U(:,1:d);
    Y=Ud'*X-Ud'*mu;
    sigma = sum(S(:,1:d));
end

