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

for i = 5  % process each subject and get feature vector for each
    t = tlist(i);
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    mammoimgright = flipdim(mammoimgright,2);

    % process images
    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.08,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.08,0);
    % get feature vector - make sure later to only determine one breast
    % cancer
    [featR] = extract_feat(processR,pecR,breastmaskR,4,15)
    [featL] = extract_feat(processL,pecL,breastmaskL,4,15)

end

