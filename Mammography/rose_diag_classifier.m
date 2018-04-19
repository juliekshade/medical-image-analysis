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

for i =1 % use 17 subjects for training and 4 for testing later (outer CV loop) 
    % just read in all training images one at a time- apply and save gabor
    % output?
    t = train(i)
    mammoimgleft(:,:) = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright(:,:) = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    % gabor filter on one at a time
    for k = orientations(1)
        j = 1
        for scale = S
            [mag(:,:,j), phase(:,:,j)] = imgaborfilt(mammoimgleft,scale,k);
            j = j+1
        end
        
    end
    %[magR, phaseR] = imgaborfilt(mammoimgright,g(1));
end


% test filter on 1st image only