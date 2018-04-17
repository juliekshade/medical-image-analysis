clear all
close all

ima = imread('img_a.gif');
imb = imread('img_b.gif');

% define training data (3 regions per class), spearately b/c confusing
csftrain = zeros(size(ima));
gmtrain = zeros(size(ima));
wmtrain = zeros(size(imb));

for i = 1:3
    csftrain = csftrain + roipoly(ima);
end

for i = 1:3
    gmtrain = gmtrain + roipoly(ima);
end

for i = 1:3
    wmtrain = wmtrain + roipoly(ima);
end

% Create an RGB image from im
imRGB(:,:,1) = ima;
imRGB(:,:,2) = ima;
imRGB(:,:,3) = ima;
greenChannel = ima;
redChannel = ima;
blueChannel = ima;
greenChannel(find(csftrain>0)) = 255;
redChannel(find(gmtrain>0)) = 255;
blueChannel(find(wmtrain>0)) = 255;
imRGB(:,:,1) = redChannel;
imRGB(:,:,3) = blueChannel;
imRGB(:,:,2) = greenChannel;
imshow(imRGB)

% plot histograms of training data for each class
imcsf = ima(find(csftrain>0));
imgm = ima(find(gmtrain>0));
imwm = ima(find(wmtrain>0));
figure(2)
histogram(imcsf(imcsf~=0),[0:1:255])
title('Histogram of CSF Training Data');
xlabel('Intensity');
ylabel('Count');
figure(3)
histogram(imgm(imgm~=0),[0:1:255])
title('Histogram of GM Training Data');
xlabel('Intensity');
ylabel('Count');
figure(4)
histogram(imwm(imwm~=0),[0:1:255])
title('Histogram of WM Training Data');
xlabel('Intensity');
ylabel('Count');

% perform k-NN classification
% aggregate training data
training = [imcsf(imcsf~=0), ones(size(imcsf(imcsf~=0)));...
    imgm(imgm~=0), ones(size(imgm(imgm~=0))).*2;
    imwm(imwm~=0), ones(size(imwm(imwm~=0))).*3]

% for each point in image b, find k nearest neighbors in training data
class = zeros(256,256,size(1:4:25,2));
kvec = 1:4:25
for k = 1:4:25
    for i = 1:256
        for j = 1:256
            if imb(i,j)~=0
                knn = (imb(i,j)-training(:,1)).^2;
                [sorted,ind] = sort(knn,'ascend');
                knn_classes = training(ind(1:k),2);
                class(i,j,(k+3)/4) = mode(knn_classes);
            end
        end
    end
end


figure(5)
imshow(class(:,:,1),[0 3])
title('Image classified with 100-NN')

% supervised gaussian classification (density est, post prob, classify)
% estimate parameters for each class from training data
u_csf = mean(imcsf(imcsf~=0))
u_gm = mean(imgm(imgm~=0))
u_wm = mean(imwm(imwm~=0))

sigma_csf = var(double(imcsf(imcsf~=0))) % dim = 1x1 because only 1 feature = intensity
sigma_gm = var(double(imgm(imgm~=0)))
sigma_wm = var(double(imwm(imwm~=0)))

P = [histcounts(training(:,2), [1 2 3 4])./size(training,1)]
P_csf = P(1)
P_gm = P(2)
P_wm = P(3)

% calculate posterior probabilties for each class for each pixel
post_prob = zeros(256,256,3)
for i = 1:256
    for j = 1:256
        if imb(i,j)~=0
        d = normpdf(double(imb(i,j)),u_csf,sigma_csf)*P_csf + ...
            normpdf(double(imb(i,j)),u_gm,sigma_gm)*P_gm + ...
            normpdf(double(imb(i,j)),u_wm,sigma_wm)*P_wm;
        post_prob(i,j,1) = (normpdf(double(imb(i,j)),u_csf,sigma_csf)*P_csf)/d;
        post_prob(i,j,2) = (normpdf(double(imb(i,j)),u_gm,sigma_gm)*P_gm)/d;
        post_prob(i,j,3) = (normpdf(double(imb(i,j)),u_wm,sigma_wm)*P_wm)/d;
        end
    end
end

% density estimate for training data
x = 1:1:255;
figure(6) 
plot(x,normpdf(x,u_csf,sigma_csf))
hold on
plot(x,normpdf(x,u_gm,sigma_gm))
plot(x,normpdf(x,u_wm,sigma_wm))
xlim = [0 255];
title('Gaussian Density Estimate for Training Data')
xlabel('Intensity')
ylabel('P(Intensity|Class)')
legend('CSF','GM','WM')

% posterior probabilites for each class
figure(7)
imshow(post_prob(:,:,1),[0 max(max(post_prob(:,:,1)))])
title('P(CSF|intensity)')

figure(8)
imshow(post_prob(:,:,2),[0 max(max(post_prob(:,:,2)))])
title('P(GM|intensity)')

figure(9)
imshow(post_prob(:,:,3),[0 max(max(post_prob(:,:,3)))])
title('P(WM|intensity)')

% find final hard classification and plot
finalclass = zeros(256,256);
for i = 1:256
    for j = 1:256
        if imb(i,j)~=0
        [m,finalclass(i,j)]=(max(post_prob(i,j,:)));
        end
    end
end

figure(10)
imshow(finalclass,[0 3]);
title('Final Gaussian Classification')





