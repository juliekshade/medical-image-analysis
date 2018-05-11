
function [mask,diagL,diagR] = mass_seg(processL,processR,pecR,pecL,breastmaskL,breastmaskR)
    
    % pre-process images by masking out pectoral muscle
    pecmaskL = poly2mask([.5 pecL(1,1) .5],[.5 .5 pecL(2,2)],size(breastmaskL,1),size(breastmaskL,2));
    processL_pec = processL.*imcomplement(pecmaskL);
    processL_norm = processL_pec./max(max(processL_pec));
    pecmaskR = poly2mask([.5 pecR(1,1) .5],[.5 .5 pecR(2,2)],size(breastmaskR,1),size(breastmaskR,2));
    processR_pec = processR.*imcomplement(pecmaskR);
    processR_norm = processR_pec./max(max(processR_pec));

    % calculate histograms for each image 
    hR =  histogram(processR_norm(find(processR_norm>0)),'Normalization','pdf');
    hold on
    hL = histogram(processL_norm(find(processL_norm>0)),'Normalization','pdf');
    Lmax = max(hL.BinCounts(:,find(hL.BinEdges(1:end-1)>.25)));
    Rmax = max(hR.BinCounts(:,find(hR.BinEdges(1:end-1)>.25)));
    
    % use peak of histograms to fistinguish more suspicious mammogram
    if Lmax<Rmax
        diagL = 1;
        diagR = 0;
        mammoimg = processL;
        pecmask=pecmaskL;
        pec=pecL;
    else
        diagL = 0;
        diagR = 1;
        mammoimg = processR;
        pecmask=pecmaskR;
        pec=pecR;
    end
    mask = zeros(size(mammoimg));

    imscale=.08;
    % preprocess image for fractal dimension analysis
    mammoseg = mammoimg.*imcomplement(pecmask);
    mammoseg = mammoseg.^.8;
    [s,I] = sort(mammoseg(find(mammoseg>0)));
    mammoseg = mammoseg./s(end-100);
    mammoseg(1:pec(2,2)*3/4,:)=0;
    k=1;
    deltaN=[];
    Ival = 0:.1:1;
    for I = Ival
        [n, r] = boxcount(mammoseg.*(mammoseg>=I));
        deltaN(k) = n(1);
        k=k+1;
    end
%    plot(-gradient(deltaN));
%    plot(gradient(-gradient(deltaN)));
    [~,I1] = max(-gradient(deltaN));
    [~,I2] = min(gradient(-gradient(deltaN)));
    lesionmask = mammoseg>(I2*.1);
    if I2<10
        se=strel('disk',5);
        CC = bwconncomp(imerode(lesionmask,se),4);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idxlist] = sort(numPixels,'descend');
        if size(numPixels,2) == 0
            diagL = 0;
            diagR = 0;
            mask = zeros(size(mammoimg));
            return 
        else
            idx = 1;
            while size(intersect(find(imdilate(pecmask,strel('disk',6))==1),...
                    CC.PixelIdxList{idxlist(idx)}),1) ~= 0
                idx = idx+1;
                if idx > size(idxlist)
                    diagL = 0;
                    diagR = 0;
                    mask = zeros(size(mammoimg));
                    return 
                end
            end
            lesionmask_biggestcc = zeros(size(lesionmask));
            lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;

            stats = regionprops(CC,'BoundingBox','Area',...
                'MajorAxisLength','MinorAxisLength');
            areabb = stats(idxlist(idx)).MajorAxisLength *...
                stats(idxlist(idx)).MinorAxisLength;
            areacc = stats(idxlist(idx)).Area;
            Aratio = areacc/areabb;
            while Aratio < .55
                idx = idx+1;
                if idx > size(idxlist)
                    idx = idx - 1; % if no cc has this property, use highest 
                end
                lesionmask_biggestcc = zeros(size(lesionmask));
                lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;
                areabb = stats(idxlist(idx)).MajorAxisLength *...
                    stats(idxlist(idx)).MinorAxisLength;
                areacc = stats(idxlist(idx)).Area;
                Aratio = areacc/areabb;
                if idx > size(idxlist)
                    Aratio = 1; 
                end
            end
        x1 = int16(stats(idxlist(idx)).BoundingBox(1));
        x2 = int16(stats(idxlist(idx)).BoundingBox(1)+...
            stats(idxlist(idx)).BoundingBox(3));
        y1 = int16(stats(idxlist(idx)).BoundingBox(2));
        y2 = int16(stats(idxlist(idx)).BoundingBox(2)+...
            stats(idxlist(idx)).BoundingBox(4));
        tumoredgemask = zeros(size(mammoseg));
        t = 30; % 30 w/ L good for tight seg

        mammoimg2 = adapthisteq(mammoimg);
        [L,NumLabels] = superpixels(mammoimg2(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
            max(x1-t,1):min(x2+t,size(mammoseg,2))),10, 'Compactness',1);
        tumoredgemask(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
            max(x1-t,1):min(x2+t,size(mammoseg,2))) = boundarymask(L,4);
        tumoredgemask = imfill(tumoredgemask,'holes');
        tumoredgemask = imopen(tumoredgemask, strel('disk',10));
        tumoredgemask = tumoredgemask.* imdilate(lesionmask_biggestcc,strel('disk',t/2));
        tumoredgemask = imdilate(tumoredgemask,strel('disk',...
            round(sum(sum(tumoredgemask))*.003)));
        expandedmask=tumoredgemask;
        if Lmax<Rmax
            mask = imresize(padarray(expandedmask,[8 8],0,'both'),1/.08);
        else
            mask = imresize(padarray(expandedmask,[8 8],0,'both'),1/.08);
            mask = flipdim(mask,2);
        end
        t = 250; % changed from 30 for feature optimization
        k=1;
        deltaN2=[];
        Ival = 0:.1:1;
        for I = Ival
            [n, r] = boxcount(mammoimg(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
            max(x1-t,1):min(x2+t,size(mammoseg,2))).*...
            (mammoseg(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
            max(x1-t,1):min(x2+t,size(mammoseg,2)))>=I));
            deltaN2(k) = n(1);
            feat2(k) = n(2)/n(1);
            k=k+1;
        end

        w_all = [61.2841   -1.8447 -236.6318];
        Xt = [1 1-deltaN2(6)/deltaN2(1) feat2(4)];

        if Xt*w_all' > 0.4848
            diagR = diagR*2;
            diagL = diagL*2;
        end
    end
end
    



%function mask = mass_seg(mammoimg,pec,breastmask)


