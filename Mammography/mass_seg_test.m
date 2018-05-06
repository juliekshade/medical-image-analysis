% find most likely thing to be a mass in each given image.
clear all;
clc;
warning('off', 'all')

datatopdir = './MammoTraining/';  
sublistfile = fullfile(['./Project1List.xlsx']);

[~,~,alllist] = xlsread(sublistfile);
sublist = alllist(2:end,1);
sublist = num2str(cell2mat(sublist));
numsubs = length(sublist);
truediag = alllist(2:end,2:3);
truediag = cell2mat(truediag);
rng(1);

for i = 8:21 % process each subject and get feature vector for each
    %  part below will be passed to function as imput
    mammoimgleft = imread([datatopdir,sublist(i,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(i,:) '_RIGHT.png']);
    mammoimgright = flipdim(mammoimgright,2);
    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.08,1);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.08,1);
    % input is thing above

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
     if truediag(i,1)>truediag(i,2)
        mammoimg = processL;
        truemaskfile = [datatopdir sublist(i,:) '_LEFT_MASK.png'];
        truemask = imread(truemaskfile);
        pecmask=pecmaskL;
        pec=pecL;
    else
        mammoimg = processR;
        truemaskfile = [datatopdir sublist(i,:) '_RIGHT_MASK.png'];
        truemask = imread(truemaskfile);
        truemask = flipdim(truemask,2);
        pecmask=pecmaskR;
        pec=pecR;
    end
     imscale=.08;
    truemask = double(imresize(truemask, imscale));
   truemask_scale = mat2gray(truemask(imscale*100:end-imscale*100,...
        imscale*100:end-imscale*100)); 
    figure(i)
    B = labeloverlay(adapthisteq(mammoimg,'Distribution','exponential'),truemask_scale);
    imshow(B)
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
    else
        idx = 1
    while size(intersect(find(imdilate(pecmask,strel('disk',6))==1),CC.PixelIdxList{idxlist(idx)})) ~= 0
        idx = idx+1
    end
    lesionmask_biggestcc = zeros(size(lesionmask));
    lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;
%     se=strel('disk',5)
%     lesionmaskfinal = imdilate(lesionmask_biggestcc,se)
    figure(i)
    imshow(imoverlay(B,lesionmask_biggestcc));
    stats = regionprops(CC,'BoundingBox','Area','MajorAxisLength','MinorAxisLength');
        areabb = stats(idxlist(idx)).MajorAxisLength * stats(idxlist(idx)).MinorAxisLength;
    areacc = stats(idxlist(idx)).Area;
    Aratio(i) = areacc/areabb;
    while Aratio(i) < .55
        idx = idx+1
        Aratio(i)
            lesionmask_biggestcc = zeros(size(lesionmask));
    lesionmask_biggestcc(CC.PixelIdxList{idxlist(idx)})=1;
    B2 = imoverlay(B,imdilate(lesionmask_biggestcc,se));
    imshow(B2)
    areabb = stats(idxlist(idx)).MajorAxisLength * stats(idxlist(idx)).MinorAxisLength;
    areacc = stats(idxlist(idx)).Area;
    Aratio(i) = areacc/areabb;
    end
    
    end
    end
%     % use lesionmask_biggest cc for classification..??
%     lesionmean = histogram(mammoseg(find(lesionmask_biggestcc.*mammoseg)>0));
%     lesionmean = I2/10
%     thresh=.01;
%     dillesion = imdilate(lesionmask_biggestcc,strel('disk',15));
%     expandedmask = zeros(size(mammoseg));
%   %  expandedmask(dillesion.*(mammoseg>(lesionmean-thresh))>0) = 1;
%   %  expandedmask(dillesion.*(mammoseg>(lesionmean-thresh))>0)=1
%     imshow(imoverlay(B,expandedmask))
   % find thresh to maximize dice coeff!
   imshow(edge(mammoseg))
   x1 = int16(stats(idxlist(idx)).BoundingBox(1))
   x2 = int16(stats(idxlist(idx)).BoundingBox(1)+stats(idxlist(idx)).BoundingBox(3))
   y1 = int16(stats(idxlist(idx)).BoundingBox(2))
   y2 = int16(stats(idxlist(idx)).BoundingBox(2)+stats(idxlist(idx)).BoundingBox(4))
   tumoredgemask = zeros(size(mammoseg))
   t = 30 % 30 w/ L good for tight seg
%  tumoredgemask(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
 %      max(x1-t,1):min(x2+t,size(mammoseg,2))) = boundarymask(mammoseg2(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
  %     max(x1-t,1):min(x2+t,size(mammoseg,2)))) 
%   imshow(imoverlay(B, tumoredgemask))
mammoimg = adapthisteq(mammoimg)
   [L,NumLabels] = superpixels(mammoimg(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
       max(x1-t,1):min(x2+t,size(mammoseg,2))),15, 'Compactness',1)
   tumoredgemask(max(y1-t,1):min(y2+t,size(mammoseg,1)),...
       max(x1-t,1):min(x2+t,size(mammoseg,2))) = boundarymask(L,4)
   tumoredgemask = imfill(tumoredgemask,'holes')
   tumoredgemask = imopen(tumoredgemask, strel('disk',3))
   tumoredgemask = tumoredgemask.* imdilate(lesionmask_biggestcc,strel('disk',t/2))
      imshow(imoverlay(B, tumoredgemask))

%   imshow(boundarymask(mammoseg(y1-10:y2+10,x1-10:x2+10)>I2/10+thresh))
end
    



%function mask = mass_seg(mammoimg,pec,breastmask)


