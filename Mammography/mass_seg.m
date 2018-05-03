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

tlist = randperm(numsubs); % random leave one out cross val ... later
for i = 7:21  % process each subject and get feature vector for each
    % boobs w only tumors!
    %  part below will be passed to function as imput
    mammoimgleft = imread([datatopdir,sublist(i,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(i,:) '_RIGHT.png']);
        mammoimgright = flipdim(mammoimgright,2);

    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.05,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.05,0);
    % input is thing above
    
    truemaskfile = [datatopdir sublist(i,:) '_LEFT_MASK.png'];
    imscale=.05
        truemask = imread(truemaskfile);
        truemask = double(imresize(truemask, imscale));
   truemask_scale = mat2gray(truemask(imscale*100:end-imscale*100,...
        imscale*100:end-imscale*100)); 
    L = 
    imshow(truemask)
    B = labeloverlay(Lenhance,truemask_scale)
    Lenhance = adapthisteq(processL);
    Renhance = adapthisteq(processR);
    figure(i+21)
    imshow(Lenhance)
    % find most likely thing to be a mass in each image? 
end
    
function mask = mass_seg(mammoimg,pec,breastmask)

    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.05,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.05,0);


mask=0

end
