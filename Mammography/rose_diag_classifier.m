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
se = strel('disk',100*imscale*5)

% loop below only collects features for the patients
for i = 1 % use 18 subjects for training and 3 for testing later (outer CV loop) 
    t = train(i)
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    % TODO: mirror image for right and make this a loop
    mammoimgleft_scale = double(imresize(mammoimgleft, imscale));
    [r,c] = size(mammoimgleft_scale)
    % try to outline breast boundary: modified traditional active contour
    g = log(1+mammoimgleft_scale);
    g_norm = mat2gray(g);
    level = graythresh(g_norm(find(g_norm>0))) % eliminate deidentification artifacts
    BW = imbinarize(g_norm,level);
    BW = imopen(BW,se); % remove text artifacts
        
    figure(i)
    imshow(g_norm, [0 1])
    hold on
    [C,h] = imcontour(BW,1,'r');
    
        % TODO: expand g_open by normal line seg analysis
    for j = 1:3
        1
    end
    N1 = [1 1]
    N2 = [1 max(find(h.ZData(1,:)>0))]
    N5 = [max(find(h.ZData(:,1)>0)) 1]
    N3 = [N5(1)/2 1]
    N4 = [N3(1) N2(2)]
    % find hough transform
    pecROI = g_norm(1:N4(1),1:N4(2));
    [H,T,R] = hough(pecROI,'Theta',30:1:80);
    P  = houghpeaks(H,1);

lines = houghlines(pecROI,T,R,P,'FillGap',5,'MinLength',sqrt(2)*N4(2));
%figure(3), imshow(pecROI,[0 1]), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
    % crop image to ROI

    % find breast mask and crop breast ?
  %  max(h.Ydata(:))

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
