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

for i = 1:21  % process each subject and get feature vector for each
    mammoimgleft = imread([datatopdir,sublist(i,:) '_LEFT.png']);
    mammoimgright = imread([datatopdir,sublist(i,:) '_RIGHT.png']);
    mammoimgright = flipdim(mammoimgright,2);

    % process images
    
    [processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.1,0);
    [processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.1,0);

    % get feature vector - make sure later to only determine one breast
    % cancer
    tic
    featR3_9_10_16(i,:) = extract_feat(processR,pecR,breastmaskR,3,20)
    featL3_9_10_16(i,:) = extract_feat(processL,pecL,breastmaskL,3,20)
    elapsedtime3_9_10_16(i) = toc

end

% step1: is one of the boobs cancer boob?
% dis part for training
% are the boobs different or not?

featL = featL3_9_10_16;
featR = featR3_9_10_16;
% find median of hist (sum up to >=.5) 
for j = 1:21
    s=0;
    k=1;
    while s<.5
        s=s+featR(j,k);
        k=k+1;
    end
    medianR(j)=k-1;
end
for j = 1:21
    s=0;
    k=1;
    while s<.5
        s=s+featL(j,k);
        k=k+1;
    end
    medianL(j)=k-1;
end
X(1,:) = abs(medianL-medianR); % X1 = difference between median angle
X(2,:) = abs(var(featR,0,2) - var(featL,0,2)); % X2 = diff between var of dist
X(3,:) = sum((featL-featR).^2,2) ;% X3 = SSD
X(4,:) = max((featL-featR).^2,[],2); % X4 = maximum squared difference
[mL,IL] = max(featL,[],2);
[mL,IR] = max(featR,[],2);
X(5,:) = abs(IL-IR); % X5 = difference between maximums of each dist
X(6,:) = [abs(median(featL,2)-median(featR,2))]'; % difference between medians
X(7,:) = abs(skewness(featL,[],2)-skewness(featR,[],2))'; % difference between skewness

X_norm = (X - min(X,[],2))./(max(X,[],2) - min(X,[],2)); % normalize feature vector

% test which features are most significant 
cancerdiag = (sum(truediag,2) >= 1);

for i = 1:1000
    t = randperm(21);
    Mdl = fitctree(X(:,t(1:18))',cancerdiag(t(1:18)));
    imp(i,:) = predictorImportance(Mdl);
    err(i) = sum(abs(predict(Mdl,X(:,t(19:21))') - cancerdiag(t(19:21)))>0)/3;
end

sum(imp>0)/1000
mean(imp)
mean(err)