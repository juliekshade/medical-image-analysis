
clear all
clc;
datatopdir = './MammoTraining/';  
sublistfile = fullfile(['./Project1List.xlsx']);
rng(1);

[~,~,alllist] = xlsread(sublistfile);
sublist = alllist(2:end,1);
sublist = num2str(cell2mat(sublist));
numsubs = length(sublist);
truediag = alllist(2:end,2:3);
truediag = cell2mat(truediag);

train = randperm(numsubs);
imscale = .15;
se = strel('disk',imscale*260);
interplinelen = imscale*1200;
bins = 0:.05:1;
tic;
% loop below only collects features for the patients
for i = 21 % use 18 subjects for training and 3 for testing later (outer CV loop) 
    t = train(i);
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    % TODO: mirror image for right and make this a loop
    % TODO: move preprocessing stuff to another file
    
    mammoimgleft_scale = double(imresize(mammoimgleft, imscale));
    mammoimgleft_scale = mat2gray(mammoimgleft_scale(imscale*100:end-imscale*100,...
        imscale*100:end-imscale*100)); % crop edges b/c artifacts 
    [r,c] = size(mammoimgleft_scale);
    % try to outline breast boundary: modified traditional active contour
    g = log(1+mammoimgleft_scale);
    g_norm = mat2gray(g);
    level = graythresh(g_norm(find(g_norm>0))); % eliminate deidentification artifacts
    BW = imbinarize(g_norm,level);
    BW = imopen(BW,se); % remove text artifacts
    stats = regionprops(BW,'Extrema'); % make sure only one region contoured or join 2 if mult
    stats = cat(2, stats.Extrema);
   	if size(stats,2) > 2  % segmented into multiple regions - start w 2, assume biggest 2
        pt1 = ceil(stats(5,1:2));
        pt2 = ceil(stats(3,3:4));
        pt3 = ceil(stats(8,3:4));
        BW = BW + poly2mask([.5 pt1(1) pt2(1) pt2(1) .5],...
        [pt1(2)-1 pt1(2)-1 pt2(2) pt3(2) pt3(2)], r, c)~=0;
    end

    figure(i) % show boob to check
    imshow(mammoimgleft_scale, [0 1])
    hold on
    [C,h] = imcontour(BW,1,'r');
    
    cx = C(1,2:end);
    cy = C(2,2:end);
    npts = imscale*500;
    ptspace = int16(C(2,1)/npts);
    boundapprox = zeros(npts,2);
    for j = 1:ptspace:C(2,1)-3 % normal line segment analysis
        x1 = cx(j);
        x2 = cx(j+2);
        y1 = cy(j);
        y2 = cy(j+2);
        mid = [(x1+x2)/2, (y1+y2)/2];
        len = pdist([x1,y1;x2,y2],'euclidean')/2;
        h = sqrt(interplinelen^2 + len^2);
        if x2-x1 >= 0 & y2-y1 >= 0
            y3 = mid(2) + interplinelen*cos(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
            x3 = mid(1) + interplinelen*sin(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 < 0 & y2-y1 >= 0
            y3 = mid(2) - interplinelen*cos(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
            x3 = mid(1) + interplinelen*sin(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 >= 0 & y2-y1 < 0
            y3 = mid(2) - interplinelen*cos(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
            x3 = mid(1) + interplinelen*sin(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 < 0 & y2-y1 < 0
            y3 = mid(2) - interplinelen*cos(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
            x3 = mid(1) - interplinelen*sin(-2*acos(len/h)-acos(abs(x1-x2)/(2*len)));
        end

        l = improfile(g_norm, [mid(1) x3], [mid(2) y3]);
        [m,index] = max(histcounts(l,bins));

        newinterplinelen = min(find(l<bins(index+1)));
        if sum(l<bins(index+1))==0
            newinterplinelen = max(find(l>0));
        end

        hnew = sqrt(newinterplinelen^2 + len^2);
        if x2-x1 >= 0 & y2-y1 >= 0
            boundapprox((j+ptspace-1)/ptspace,2) = mid(2) + newinterplinelen*cos(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
            boundapprox((j+ptspace-1)/ptspace,1) = mid(1) + newinterplinelen*sin(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 < 0 & y2-y1 >= 0
            boundapprox((j+ptspace-1)/ptspace,2) = mid(2) - newinterplinelen*cos(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
            boundapprox((j+ptspace-1)/ptspace,1) = mid(1) + newinterplinelen*sin(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 >= 0 & y2-y1 < 0
            boundapprox((j+ptspace-1)/ptspace,2) = mid(2) - newinterplinelen*cos(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
            boundapprox((j+ptspace-1)/ptspace,1) = mid(1) + newinterplinelen*sin(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
        elseif x2-x1 < 0 & y2-y1 < 0
            boundapprox((j+ptspace-1)/ptspace,2) = mid(2) - newinterplinelen*cos(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
            boundapprox((j+ptspace-1)/ptspace,1) = mid(1) - newinterplinelen*sin(-2*acos(len/hnew)-acos(abs(x1-x2)/(2*len)));
        end

    end
    boundapprox = int16(boundapprox);
    boundapprox(boundapprox(:,1)<1,2)=1;
    boundapprox(boundapprox(:,1)>c,1)=c;
    boundapprox(boundapprox(:,2)<1,2)=1;
    boundapprox(boundapprox(:,2)>r,1)=r;
    if boundapprox(1,2) ~= 1
        boundapprox(1,2) = 1;
        boundapprox(1,1) = boundapprox(find(min(boundapprox(:,2))),1);
    end
    [C,ia,ic] = unique(boundapprox(:,1),'stable');
    boundapprox = double(boundapprox(ia,:));
    boundapprox(end,1) = 1;
    boundapprox(end,2) = boundapprox(end-1,2);
    
    % add result of normal line segment analysis to plot
    scatter(boundapprox(:,1),boundapprox(:,2));

    % create smooth contour based on normal line segment analysis
    windowWidth = imscale*160-1;
    polynomialOrder = 2;
    smoothX = [sgolayfilt(boundapprox(:,1), polynomialOrder, windowWidth); 1];
    smoothY = [1; sgolayfilt(boundapprox(:,2), polynomialOrder, windowWidth)];
    smoothX = [smoothX(1); smoothX];
    smoothY = [smoothY; smoothY(end)];
    
    % plot this smooth contour
    plot(smoothX,smoothY);
    
    % determine critical points for hough transform
    N1 = [1 1];
    N2 = [1 smoothX(1,1)];
    N5 = [smoothY(end) 1];
    N3 = [N5(1)*2/3 1];
    N4 = [N3(1) N2(2)];
    
    % find hough transform
    pecROI = g_norm(1:N4(1),1:N4(2));
    [H,T,R] = hough(pecROI,'Theta',20:1:80);
    P  = houghpeaks(H,100);
    lines = houghlines(pecROI,T,R,P);
    % remove lines that don't intersect N2-N3
    toremove = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       if xy(1,2) ~= 1 || xy(2,1) ~=1
           toremove = [toremove k];
       end
    end
    lines(toremove) = [];
    bestline=1;
    intensitymax = 0;
    for k = 1:length(lines)     
        xy = [lines(k).point1; lines(k).point2];
        trimask = poly2mask([.5 xy(:,1)'],[.5 xy(:,2)'],r,c);
        intensity = mean(g_norm(trimask==1))/var(g_norm(trimask==1));
        if intensity > intensitymax
            bestline=k;
            intensitymax = intensity;
        end
    end
    xy = [lines(bestline).point1; lines(bestline).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
    % Plot beginnings and end of best line
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
    % find breastMask
    breastMask = zeros(r,c);
    breastMask(poly2mask([.5; smoothX],[.5; smoothY],r,c))=1;
    
    maskBreast = g_norm.*breastMask;
   %TODO: density correction in breast margin
    D = bwdist(imcomplement(breastMask))
    imshow(D,[0 90])
    hold on
    [C,h] = imcontour(D,10)
end
toc;
% RETURN: breast mask, corrected image, and pectoral muscle line?


%     % gabor filter on one at a time
%     for k = orientations(1)
%         j = 1
%         for scale = S
%             [mag(:,:,j), phase(:,:,j)] = imgaborfilt(mammoimgleft,scale,k);
%             j = j+1
%         end
%         
%     end
%     %[magR, phaseR] = imgaborfilt(mammoimgright,g(1));


