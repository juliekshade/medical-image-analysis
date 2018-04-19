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

    t = train(i)
    mammoimgleft(:,:) = imread([datatopdir,sublist(t,:) '_LEFT.png']);
    mammoimgright(:,:) = imread([datatopdir,sublist(t,:) '_RIGHT.png']);
    % segment breast region first? or crop and enhance
    % try hough transform
    % assume at least 10x10 px of top left are pectoral
    x = size(mammoimgleft
    N1 = [0 0]
    N2 = [0 min(find(mammoimgleft(50,:)<mean(mean(mammoimgleft(50:100,50:100)))/3.0))]
    N3 = [min(find(mammoimgleft(:,50)<mean(mean(mammoimgleft(50:100,50:100)))/3.0)) 0]
    
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