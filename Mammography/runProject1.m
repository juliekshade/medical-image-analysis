 function [estdiag, estmaskleft,estmaskright] = runProject1(mammoimgleft,mammoimgright)
% Inputs:   mammoimgleft  -     the mammogram of the left side (rowL x colL)
%                               where rowL and colL are the image
%                               dimensions.
%           mammoimgright -     the mammogram of the right side (rowR x colR)
%
% Outputs:  estdiag      -   the estimated diagnosis.  Should only have
%                            values of 0(healthy), 1(benigh), and 2(cancer).
%                            and be of size (1 x 2) (left, right).
%           estmaskleft  -   the output binary mask for the left side. 
%                            Should only have values of zero or one.
%                            If the estdiag is 0, the mask should only have 
%                            values of 0. Should be of size (rowL x colL).
%           estmaskright -   the output binary mask for the right side. 
%                            Should only have values of zero or one.
%                            If the estdiag is 0, the mask should only have 
%                            values of 0. Should be of size (rowR x colR).
%

% MAKE UP SOME BAD RESULTS FOR TESTING :)
estdiag = zeros(1,2);
estmaskleft = zeros(size(mammoimgleft));
estmaskright = zeros(size(mammoimgright));
warning('off','all')

%% PUT IN YOUR DIAGNOSIS AND SEGMENTATION CODE BELOW!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My tesging code randomly diagnose the images and labels 10% of the pixels
% as 1 if the diagnosis is not healthy.
rr = size(mammoimgright,1);
rc = size(mammoimgright,2);
lr = size(mammoimgleft,1);
lc = size(mammoimgleft,2);
% pre-process images
mammoimgright = flipdim(mammoimgright,2);
[processR,pecR,breastmaskR] = mammo_preprocess(mammoimgright,.08,0);
[processL,pecL,breastmaskL] = mammo_preprocess(mammoimgleft,.08,0);

% classify patient as cancer or no cancer?
iscancer = classify_cancer(processL,processR,pecR,pecL,breastmaskL,breastmaskR);

% if cancer, try to segment mass and return diagnosis
if iscancer==1
    [mask,diagL,diagR] = mass_seg(processL,processR,pecR,pecL,breastmaskL,breastmaskR);
    if diagL >0
        mask = imresize(mask,[lr lc]);
    else
        mask = imresize(mask,[rr rc]);
    end
    estdiag = [diagL diagR]
else
    estdiag = [0 0]
end

% Left side
if estdiag(1) ~= 0
    % TODO: find estimated mask for left side breast
    estmaskleft = mask;
end
% Right side
if estdiag(2) ~= 0
    % TODO: find estimated mask for right side breast
    estmaskright = mask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



