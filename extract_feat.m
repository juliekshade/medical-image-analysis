% feature extraction - histogram of angular information
function [feat] = extract_feat(mammoimg, pec, breastmask, numscales, dtheta)
    % crop closer to breast region
    mammoimg_crop = mammoimg(pec(2,2)/2:max(sum(breastmask,1)),1:max(sum(breastmask,2)));
    breastmask_crop = breastmask(pec(2,2)/2:max(sum(breastmask,1)),1:max(sum(breastmask,2)));

    r = size(mammoimg_crop,1);
    c = size(mammoimg_crop,2);
    wavelengthMin = 4/sqrt(2);
    wavelengthMax = hypot(r,c);
    n = floor(log2(wavelengthMax/wavelengthMin));
    wavelength = 1.6.^(0:n/(numscales-1):n) * wavelengthMin;

    orientation = 0:dtheta:(180-dtheta);

    g = gabor(wavelength,orientation);
    gabormag = imgaborfilt(mammoimg_crop,g);

    k_imgs = zeros(r*c,size(orientation,2));
    for k = 1:size(orientation,2)
        a = reshape(gabormag(:,:,(k*numscales)-numscales+1:k*numscales),[r*c,numscales]);
        [coeff,score,~,~,explained,~] = pca(a);
        sum_var = 0;
        it=0;
        while sum_var < 99
            sum_var = sum_var+explained(it+1);
            it = it+1;
        end
        Xcentered = score*coeff';
        k_imgs(:,k) = sum(Xcentered(:,1:it),2);
    end
    
    % binarize images in mask area only - elim edge artifacts by mask
    D = bwdist(imcomplement(breastmask_crop));
    breastmask_crop = breastmask_crop.*(D>6);
    kthresh = zeros(r,c,size(orientation,2));
    for i = 1:size(orientation,2)
        thisimg = reshape(k_imgs(:,i),[r c]).*breastmask_crop;
        [thresh(i),EM(i)] = graythresh(thisimg);
        k_thresh(:,:,i) = thisimg;
    end
    [m,I] = max(EM);
    k_thresh=k_thresh>thresh(I);
%     figure;
%     imshow(sum(k_thresh,3));
    s = sum(sum(k_thresh,1),2)./sum(sum(sum(k_thresh,1),2),3);
    feat = s(:);
end
