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

train = randperm(numsubs)
orientations = [0:180/11:180]
S = [4:2:10]
imscale = .2
se = strel('disk',50*imscale*5)
interplinelen = imscale*1000
bins1 = 0:.01:1
bins = 0:.05:1

% loop below only collects features for the patients
for i = 8 % use 18 subjects for training and 3 for testing later (outer CV loop) 
    t = train(i)
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    % TODO: mirror image for right and make this a loop
    % TODO: move preprocessing stuff to another file
    mammoimgleft_scale = double(imresize(mammoimgleft, imscale));
    mammoimgleft_scale = mat2gray(mammoimgleft_scale(imscale*100:end-imscale*100,imscale*100:end-imscale*100)) % crop edges b/c artifacts 
    [r,c] = size(mammoimgleft_scale)
    % try to outline breast boundary: modified traditional active contour
    g = log(1+mammoimgleft_scale);
    g_norm = mat2gray(g);
    level = graythresh(g_norm(find(g_norm>0))) % eliminate deidentification artifacts
    BW = imbinarize(g_norm,level);
    BW = imopen(BW,se); % remove text artifacts
        
    figure(i)
    imshow(mammoimgleft_scale, [0 1])
    hold on
    [C,h] = imcontour(BW,1,'r');
    
    cx = C(1,2:end);
    cy = C(2,2:end);
    boundapprox = []
    npts = 30
    ptspace = int16(C(2,1)/npts)
  %  ptspace = imscale*70
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

        hnew = sqrt(newinterplinelen^2 + len^2);
        (j+ptspace-1)/ptspace
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
        boundapprox(1,2) = 1
        boundapprox(1,1) = boundapprox(2,1)+1
    end
    [C,ia,ic] = unique(boundapprox(:,1),'stable')


    boundapprox = double(boundapprox(ia,:))
    boundapprox(end,1) = 1

    scatter(boundapprox(:,1),boundapprox(:,2))

    N1 = [1 1];
    N2 = [1 boundapprox(1,1)]
    N5 = [boundapprox(end,2) 1]
    N3 = [N5(1)*1/2 1];
    N4 = [N3(1) N2(2)];
    
    % find hough transform
    pecROI = g_norm(1:N4(1),1:N4(2));
    [H,T,R] = hough(pecROI,'Theta',20:1:80);
    P  = houghpeaks(H,100);
    lines = houghlines(pecROI,T,R,P,'MinLength',double(sqrt(2)*N4(2)));
    % remove lines that don't intersect N2-N3
    toremove = []
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       if xy(1,2) ~= 1
           toremove = [toremove k]
       end
    end
    lines(toremove) = [];
    bestline=1;
    intensitymax = 0;
    for k = 1:length(lines)     
        xy = [lines(k).point1; lines(k).point2];
        trimask = poly2mask([1 xy(:,1)'],[1 xy(:,2)'],r,c);
        intensity = mean(g_norm(trimask==1))/var(g_norm(trimask==1))*sum(sum(trimask));
        if intensity > intensitymax;
            bestline=k;
            intensitymax = intensity;
        end
    end
    xy = [lines(bestline).point1; lines(bestline).point2];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
   % Plot beginnings and end of best line
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
   % crop image to breastROI
   breastROI = mammoimgleft_scale((xy(2,1)+xy(2,2))/2:max(boundapprox(:,2)), 1:max(boundapprox(:,1)))

 %  pp = spline(1:c/size(boundapprox,1):c,boundapprox');
      pp = csaps(boundapprox(:,1),boundapprox(:,2), .4);

   yy = ppval(pp, boundapprox(:,1));
 %  plot(yy(1,:),yy(2,:),'-b')
   plot(boundapprox(:,1),yy,'-b')
%    breastMask = zeros(r,c)
%    breastMask(poly2mask([1; boundapprox(:,1)],[1; yy],r,c))=1;
   breastMask(poly2mask([1 xy(:,1)'],[1 xy(:,2)'],r,c))=0;
    breastMask(poly2mask([xy(:,1)'+10],[xy(:,2)'+10],r,c))=3
   figure(i+1)
   imshow(breastROI)
   % get rid of pectroal muscle
   
%    figure(i+2)
%   % imshow(breastMask)
% 
%     breastROImask = breastMask((xy(2,1)+xy(2,2))/2:max(boundapprox(:,2)), 1:max(boundapprox(:,1)))
%     maskedBreastROI = breastROI.*breastROImask
%     imshow(maskedBreastROI)
end


%     % segment breast region first? or crop and enhance
%     % try multiple image scales?
%     % try hough transform
%     [r, c] = size(mammoimgleft)
%     re = int16(r/50)
%     ce = int16(c/50)
%     N1 = [0 0]
%     N2 = [0 min(find(mammoimgleft(re,:)<mean(mean(mammoimgleft(re:2*re,ce:ce*2)))/6))] % breast bound
%     N3 = [min(find(mammoimgleft(:,ce)<mean(mean(mammoimgleft(re:2*re,ce:ce*2)))/2)) 0] % pec muscle bound
%     
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


% here is where the actual classification stuff starts
