% training script project 1 
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
for i = 1:21  % process each subject and get feature vector for each
    %  part below will be passed to function as imput
    t = tlist(i);
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.05,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.05,0);
    % find most likely thing to be a mass in each image?
end
    
function mask = mass_seg(mammoimg,pec,breastmask)

    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.05,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.05,0);


mask=0

end
