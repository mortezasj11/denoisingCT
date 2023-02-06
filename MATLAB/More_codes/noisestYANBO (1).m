% Convolution and edge exclusion based  standard derivertive estimation
%
% 2016-10-11

function est_sd = noisestYANBO(I0)

p = 0.9;    % percentage of avaiable pixels

[sx, sy, nMC] = size(I0);
filterxx = [-1 -2 -1; 0 0 0; 1 2 1];        % Sobel along x and y
filteryy = [-1 0 1; -2 0 2; -1 0 1];

filterxx = filterxx*sqrt(6)/12;     % normallization
filteryy = filteryy*sqrt(6)/12;

% compute the difference image of the sum image
imsumxx = imfilter(sum(I0, 3), filterxx);
imsumyy = imfilter(sum(I0, 3), filteryy);
Isumdiff = imsumxx.^2 + imsumyy.^2;

% get the location of availabel pixels
% Imask1 = findthresh(reshape(Isumdiff, sx*sy, 1), p);
Isumdiffsort1 = sort(reshape(Isumdiff, sx*sy, 1));
pos = round(length(Isumdiffsort1)*p);
threshp = Isumdiffsort1(pos);
Imask = reshape(Isumdiff<threshp, sx, sy);

est_sd2 = 0;
% difference of the each channel image
for iMC = 1:nMC
    imxx = imfilter(I0(:,:,iMC), filterxx);
    imyy = imfilter(I0(:,:,iMC), filteryy);
    Idiff = imxx.^2 + imyy.^2;    
    est_sd2 = est_sd2 + sum(sum(Idiff.*Imask));    
end

est_sd = sqrt(est_sd2/(sx*sy*p)/nMC);

end

