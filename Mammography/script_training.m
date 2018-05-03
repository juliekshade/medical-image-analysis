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
    t = tlist(i);
    mammoimgleft = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    mammoimgright = flipdim(mammoimgright,2);

    % process images
    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.05,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.05,0);

    % get feature vector - make sure later to only determine one breast
    % cancer
    tic
    featR415_05(i,:) = extract_feat(processR,pecR,breastmaskR,4,15)
    featL415_05(i,:) = extract_feat(processL,pecL,breastmaskL,4,15)
    elapsedtime05(i,1) = toc
    tic
    featR420_05(i,:) = extract_feat(processR,pecR,breastmaskR,4,20)
    featL420_05(i,:) = extract_feat(processL,pecL,breastmaskL,4,20)
    elapsedtime05(i,2) = toc
    tic
    featR315_05(i,:) = extract_feat(processR,pecR,breastmaskR,3,15)
    featL315_05(i,:) = extract_feat(processL,pecL,breastmaskL,3,15)
    elapsedtime05(i,3) = toc
    tic
    featR320_05(i,:) = extract_feat(processR,pecR,breastmaskR,3,20)
    featL320_05(i,:) = extract_feat(processL,pecL,breastmaskL,3,20)
    elapsedtime05(i,4) = toc
end
for i = 1:21
    tlist(i)
    sorteddiag(i,:)=truediag(tlist(i),:) % LEFT THEN RIGHT BOOB
end
% step1: is one of the boobs cancer boob?
% dis part for training
% are the boobs different or not?

featL = featL320_05;
featR = featR320_05;
% find median of hist (sum up to >=.5) 
for j = 1:21
    sum=0
    k=1
    while sum<.5
        sum=sum+featR(j,k)
        k=k+1
    end
    medianR(j)=k-1
end
% X is feature vector - X1 = difference between medians
X(1,:) = abs(medianL-medianR)
% find stdev of hist?
var(featR,0,2)
% find SSD

for j = 1:21
    
