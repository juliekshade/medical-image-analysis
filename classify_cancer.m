% predict whether ONE of the two given breasts has a tumor (doesn't tell
% which) - note will hopefully elim false pos later
% developed using script_training for evaluation of many models

function diag = classify_cancer(processL,processR,pecR,pecL,breastmaskL,breastmaskR)

    featR = extract_feat(processR,pecR,breastmaskR,3,20);
    featL = extract_feat(processL,pecL,breastmaskL,3,20);
    
    s=0;
    k=1;
    while s<.5
        s=s+featR(k);
        k=k+1;
    end
    medianR=k-1;
    s=0;
    k=1;
    while s<.5
        s=s+featL(k);
        k=k+1;
    end
    medianL=k-1;

    w_all = [0.5024    0.3467   22.4738];
    Xt = [1 abs(medianL-medianR) sum((featL-featR).^2)];
    if Xt*w_all' > .5091
        diag = 1;
    else
        diag = 0;
    end
end