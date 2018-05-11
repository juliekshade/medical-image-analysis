function [process, pec, mask] = mammo_preprocess(mammoimg,imscale,displ)
% Inputs:   mammoimg     -     the (fullsize) mammogram to process
%           imscale      -     factor to rescale image by for processing 
%                              speed, lower values increase speed
%           displ        -     indicates whether to plot figures
%
% Outputs:  process      -   
%           pec          -   
%           mask         -   the output binary mask for the breast region. 
%                            Should only have values of zero or one.
%                            Should be of size (rowR x colR).
%

%     if displ == 0
%         set(0,'DefaultFigureVisible','off');
%     end
    se = strel('disk',round(imscale*260));
    interplinelen = imscale*1200;
    bins = 0:.05:1;

    % rescale, crop, and enhance image for processing
    mammoimg_scale = double(imresize(mammoimg, imscale));
    mammoimg_scale = mat2gray(mammoimg_scale(imscale*100:end-imscale*100,...
        imscale*100:end-imscale*100)); 
    [r,c] = size(mammoimg_scale);
    g = log(1+mammoimg_scale);
    g_norm = mat2gray(g);

    % eliminate deidentification artifacts and binarize image using Otsu
    level = graythresh(g_norm(find(g_norm>0))); 
    BW = imbinarize(g_norm,level);
    BW = imopen(BW,se); % remove text artifacts

    % sometimes selcts pectoral + breast separately: merge if this happens
    stats = regionprops(BW,'Extrema'); 
    stats = cat(2, stats.Extrema);
    if size(stats,2) > 2  
        pt1 = ceil(stats(5,1:2));
        pt2 = ceil(stats(3,3:4));
        pt3 = ceil(stats(8,3:4));
        BW = BW + poly2mask([.5 pt1(1) pt2(1) pt2(1) .5],...
        [pt1(2)-1 pt1(2)-1 pt2(2) pt3(2) pt3(2)], r, c)~=0;
    end

    % create initial contour and display to check accuracy if desired
%     h(1) = figure; % show boob to check
%     imshow(mammoimg_scale, [0 1]);
%     hold on
figure(1)
     [C,h] = imcontour(BW,1,'r');

    cx = C(1,2:end);
    cy = C(2,2:end);

    % normal line segment analysis on originial contour with imscale*500 lines
    npts = imscale*500;
    ptspace = int16(C(2,1)/npts);
    boundapprox = zeros(npts,2);
    dontinclude = [];
    for j = 1:ptspace:C(2,1)-3 
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
        if newinterplinelen < 3
            dontinclude = [dontinclude (j+ptspace-1)/ptspace];
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
    if size(dontinclude,1) > 0 
      boundapprox(dontinclude,:) = []; 
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

    [~,ia,~] = unique(boundapprox(:,1),'stable');
    boundapprox = double(boundapprox(ia,:));
    boundapprox(end,1) = 1;
    boundapprox(end,2) = boundapprox(end-1,2);

    % add result of normal line segment analysis to plot
%     scatter(boundapprox(:,1),boundapprox(:,2));

    % create smooth contour based on normal line segment analysis
    windowWidth = 2*floor(imscale*80)+1;
    polynomialOrder = 2;
    smoothX = [sgolayfilt(boundapprox(:,1), polynomialOrder, windowWidth); .5];
    smoothY = [.5; sgolayfilt(boundapprox(:,2), polynomialOrder, windowWidth)];
    smoothX = [smoothX(1); smoothX];
    smoothY = [smoothY; smoothY(end)];

    % plot this smooth contour
%     plot(smoothX,smoothY);

    % determine critical points for hough transform
    N1 = [1 1];
    N2 = [1 smoothX(1)];
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
       if xy(1,2) > 1 || xy(2,1) > 1
           toremove = [toremove k];
       end
    end
    if size(toremove) < length(lines)
        lines(toremove) = [];
    end
    
    % find best guess of line defining pectoral muscle
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
%     if displ == 1
%         plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
%         plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%         plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%     end

    % find breastMask
    breastMask = zeros(r,c);
    breastMask(poly2mask([.5; smoothX],[.5; smoothY],r,c))=1;

    % density correction in breast margin
    D = bwdist(imcomplement(breastMask));
    ncont = imscale*300;
    C = contourc(double(D), 0:1:ncont);
    contlevels = zeros(1,int16(ncont));
    distancelevels = zeros(1,int16(ncont));
    ind1 = 1;
    for j = 1:ncont
        npts = C(2,ind1);
        distancelevels(j) = C(1,ind1);
        contlevels(j) = mean(mean(mammoimg_scale(round(C(:,ind1+1:npts+ind1)))));
        ind1 = npts+ind1+1;
    end
    inside = mean(contlevels(ncont*2/3:end));
    p = polyfit(distancelevels(1:ncont),contlevels(1:ncont),8);

    % create enhancement mask
    enhancementmask = inside - polyval(p,D);
    enhancementmask(D>distancelevels(ncont)) = 0;
    mammoimgscale_enhanced = mammoimg_scale+enhancementmask;
    maskBreastEnhanced = mammoimgscale_enhanced.*breastMask;
%     h(1) = figure;
%     imshow(maskBreastEnhanced)

    mask = breastMask;
    process = maskBreastEnhanced;
    pec = xy;
    
%     if displ == 1
%         figure(1)
%         figure(2)
%     end
end




