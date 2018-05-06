% predict whether ONE of the two given breasts has a tumor (doesn't tell
% which) - note will hopefully elim false pos later
% developed using script_training for evaluation of many models

function diag = classify_cancer(processL,processR,pecR,pecL,breastmaskL,breastmaskR)

    featR = extract_feat(processR,pecR,breastmaskR,3,20);
    featL = extract_feat(processL,pecL,breastmaskL,3,20);
    
    SSD = sum((featL-featR).^2);
    if SSD > .0004213
        diag = 1;
    else
        diag = 0;
    end
end