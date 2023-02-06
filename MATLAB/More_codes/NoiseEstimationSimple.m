
function Sigmahat= NoiseEstimationSimple(Image )
[H, W]=size(Image);

M=[1 -2 1; -2 4 -2; 1 -2 1];
sumImgI=sum(sum(abs(conv2(Image, M,'same'))));

% 3 Averaging
Sigmahat=sumImgI*sqrt(0.5*pi)./(6*(W-2)*(H-2));

end