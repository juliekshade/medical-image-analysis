
function [mask,diagL,diagR] = mass_seg(processL,processR,pecR,pecL,breastmaskL,breastmaskR)
    pecmaskL = poly2mask([.5 pecL(1,1) .5],[.5 .5 pecL(2,2)],size(breastmaskL,1),size(breastmaskL,2));
    processL_pec = processL.*imcomplement(pecmaskL);
    processL_norm = processL_pec./max(max(processL_pec));
    pecmaskR = poly2mask([.5 pecR(1,1) .5],[.5 .5 pecR(2,2)],size(breastmaskR,1),size(breastmaskR,2));
    processR_pec = processR.*imcomplement(pecmaskR);
    processR_norm = processR_pec./max(max(processR_pec));

   % figure(1)
    hR =  histogram(processR_norm(find(processR_norm>0)),'Normalization','pdf');
    hold on
    hL = histogram(processL_norm(find(processL_norm>0)),'Normalization','pdf');

    Lmax = max(hL.BinCounts(:,find(hL.BinEdges(1:end-1)>.25)));
    Rmax = max(hR.BinCounts(:,find(hR.BinEdges(1:end-1)>.25)));
    
%     if Lmax<Rmax
%         mammoseg = processL;
%         truemaskfile = [datatopdir sublist(i,:) '_LEFT_MASK.png'];
%         truemask = imread(truemaskfile);
%     else
%         mammoseg = processR;
%         truemaskfile = [datatopdir sublist(i,:) '_RIGHT_MASK.png'];
%         truemask = imread(truemaskfile);
%         truemask = flipdim(truemask,2);
%     end
%      if Lmax<Rmax
%          diagL = 1;
%          diagR = 0;
%         mammoimg = processL;
% %         truemaskfile = [datatopdir sublist(i,:) '_LEFT_MASK.png'];
% %         truemask = imread(truemaskfile);
%          pecmask=pecmaskL;
%         pec=pecL;
%      else
%         diagL = 0;
%         diagR=1;
%         mammoimg = processR;
% %         truemaskfile = [datatopdir sublist(i,:) '_RIGHT_MASK.png'];
% %         truemask = imread(truemaskfile);
% %         truemask = flipdim(truemask,2);
%         pecmask=pecmaskR;
%         pec=pecR;
%      end
    if Lmax<Rmax
        diagL = 1;
         diagR = 0;
        mammoimg = processL;
%         truemaskfile = [datatopdir sublist(i,:) '_LEFT_MASK.png'];
%         truemask = imread(truemaskfile);
        pecmask=pecmaskL;
        pec=pecL;
    else
        diagL = 0;
         diagR = 1;
        mammoimg = processR;
%         truemaskfile = [datatopdir sublist(i,:) '_RIGHT_MASK.png'];
%         truemask = imread(truemaskfile);
%         truemask = flipdim(truemask,2);
        pecmask=pecmaskR;
        pec=pecR;
    end
            mask = zeros(size(mammoimg));

     imscale=.08;
%     truemask = double(imresize(truemask, imscale));
%    truemask_scale = mat2gray(truemask(imscale*100:end-imscale*100,...
%         imscale*100:end-imscale*100)); 
%     figure(i)
%     B = labeloverlay(adapthisteq(mammoimg,'Distribution','exponential'),truemask_scale);
%     imshow(B)
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
   plot(-gradient(deltaN));
   plot(gradient(-gradient(deltaN)));
    [~,I1] = max(-gradient(deltaN));
    [~,I2] = min(gradient(-gradient(deltaN)));
    lesionmask = mammoseg>(I2*.1);
    if I2>=10
        disp('no lesion');
    else
        se=strel('disk',5)
    CC = bwconncomp(imerode(lesionmask,se),4);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idxlist] = sort(numPixels,'descend');
    if size(numPixels,2) == 0
        disp('no lesion');
        diagL = 0;
        diagR = 0;
        mask = zeros(size(mammoimg))
        return 
    else
        idx = 1
    while size(intersect(find(imdilate(pecmask,strel('disk',6))==1),CC.PixelIdxList{idxlist(idx)})) ~= 0
        intersect(find(imdilate(pecmask,strel('disk',6))==1),CC.PixelIdxList{idxlist(idx)})
        idx = idx+1
        if idx > size(idxlist)
            disp('no lesion');
        diagL = 0;
        diagR = 0;
        mask = zeros(size(mammoimg));
        return 
        end
    end
    lesionmask_biggestcc = zeros(size(lesionmask));
    lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;
%     se=strel('disk',5)
%     lesionmaskfinal = imdilate(lesionmask_biggestcc,se)
%     figure(i)
%     imshow(imoverlay(B,lesionmask_biggestcc));
    stats = regionprops(CC,'BoundingBox','Area','MajorAxisLength','MinorAxisLength');
        areabb = stats(idxlist(idx)).MajorAxisLength * stats(idxlist(idx)).MinorAxisLength;
    areacc = stats(idxlist(idx)).Area;
    Aratio = areacc/areabb;
    while Aratio < .55
        idx = idx+1;
        
            lesionmask_biggestcc = zeros(size(lesionmask));
    lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;
%     B2 = imoverlay(B,imdilate(lesionmask_biggestcc,se));
%     imshow(B2)
    areabb = stats(idxlist(idx)).MajorAxisLength * stats(idxlist(idx)).MinorAxisLength;
    areacc = stats(idxlist(idx)).Area;
    Aratio = areacc/areabb;
    end
%     lesionmean = mean(mammoseg(find(lesionmask_biggestcc.*mammoseg)>0));
%     thresh=.03;
%     dillesion = imdilate(lesionmask_biggestcc,strel('disk',20));
%     expandedmask = zeros(size(mammoseg));
%     expandedmask(dillesion.*(mammoseg<(lesionmean+thresh))>0) = 1;
%     expandedmask(dillesion.*(mammoseg>(lesionmean-thresh))>0)=1;
%     % find thresh to maximize dice coeff!
x1 = int16(stats(idxlist(idx)).BoundingBox(1));
   x2 = int16(stats(idxlist(idx)).BoundingBox(1)+stats(idxlist(idx)).BoundingBox(3));
   y1 = int16(stats(idxlist(idx)).BoundingBox(2));
   y2 = int16(stats(idxlist(idx)).BoundingBox(2)+stats(idxlist(idx)).BoundingBox(4));
   tumoredgemask = zeros(size(mammoseg));
    t = 30 % 30 w/ L good for tight seg
%  tumoredgemask(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
 %      max(x1-t,1):min(x2+t,size(mammoseg,2))) = boundarymask(mammoseg2(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
  %     max(x1-t,1):min(x2+t,size(mammoseg,2)))) 
%   imshow(imoverlay(B, tumoredgemask))
    mammoimg = adapthisteq(mammoimg);
   [L,NumLabels] = superpixels(mammoimg(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
       max(x1-t,1):min(x2+t,size(mammoseg,2))),10, 'Compactness',1)
   tumoredgemask(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
       max(x1-t,1):min(x2+t,size(mammoseg,2))) = boundarymask(L,4);
   tumoredgemask = imfill(tumoredgemask,'holes');
   tumoredgemask = imopen(tumoredgemask, strel('disk',3));
   tumoredgemask = tumoredgemask.* imdilate(lesionmask_biggestcc,strel('disk',t/2));
   %   imshow(imoverlay(B, tumoredgemask));
      expandedmask=tumoredgemask;
    if Lmax<Rmax
        mask = imresize(padarray(expandedmask,[8 8],0,'both'),1/.08);
    else
                mask = imresize(padarray(expandedmask,[8 8],0,'both'),1/.08);
        mask = flipdim(mask,2);
    end
    end
        
%     lesionmean = mean(mammoseg(find(lesionmask_biggestcc.*mammoseg)>0));
%     thresh=.02;
%     dillesion = imdilate(lesionmask_biggestcc,strel('disk',20));
%     expandedmask = zeros(size(mammoseg));
%     expandedmask(dillesion.*(mammoseg<(lesionmean+thresh))>0) = 1;
%     expandedmask(dillesion.*(mammoseg>(lesionmean-thresh))>0)=1;
%     % find thresh to maximize dice coeff!
%     mask = expandedmask;
end
    



%function mask = mass_seg(mammoimg,pec,breastmask)


