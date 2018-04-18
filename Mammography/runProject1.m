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

%% PUT IN YOUR DIAGNOSIS AND SEGMENTATION CODE BELOW!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My tesging code randomly diagnose the images and labels 10% of the pixels
% as 1 if the diagnosis is not healthy.
diag = [0,1,2];
estpos = randi(length(diag),[1,2]);
estdiag = diag(estpos);

% Left side
if estdiag(1) ~= 0
    % TODO: find estimated mask for left side breast
    estmaskleft( rand(size(estmaskleft)) < 0.1 ) = 1;
end
% Right side
if estdiag(2) ~= 0
    % TODO: find estimated mask for right side breast
    estmaskright( rand(size(estmaskright)) < 0.1 ) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



