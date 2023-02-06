% load('lena.mat')
% Image1=double(lena);
% percentage=20;
% 
% Sigmaa=1:40;
% for i=1:numel(Sigmaa)
%     
%     Image=Image1 + Sigmaa(i)*randn(size(Image1));
%     Sigma_hat(i) = NoiseEstimationMorii( Image, percentage );
% end
% plot(Sigmaa,Sigma_hat,'*')
% grid on


function [SigmaAfterfitting,Sigmahat]= NoiseEstimationMorii( Image, percentage )


if nargin==1
    percentage=10;
end

% 1 SobelEdgeDetection
maskX = [-1 0 1 ; -2 0 2; -1 0 1];
maskY = [-1 -2 -1 ; 0 0 0 ; 1 2 1] ;

resX = conv2(Image, maskX);
resY = conv2(Image, maskY);

magnitude=abs(resX)+abs(resY);%magnitude1 = sqrt(resX.^2 + resY.^2);%imshow(magnitude,[])

binNumber=300;
[N,edges] = histcounts(magnitude,binNumber);
Maxx=sum(N);

q=0;
for i=1:binNumber
q=q + N(i);
    if q>(Maxx/percentage);
        break
    end
end
threshold=edges(i);
magnitude(magnitude<threshold)=0;
%test=sum(N(1:i))/sum(N);%should be close to 10%

% 2 Laplacian Convolution
[H W]=size(magnitude);
magnitude=double(magnitude);

M=[1 -2 1; -2 4 -2; 1 -2 1];
sumImgI=sum(sum(abs(conv2(magnitude, M))));

% 3 Averaging
Sigmahat=sumImgI*sqrt(0.5*pi)./(6*(W-2)*(H-2));

p =[2.3416    6.9744];
SigmaAfterfitting=(Sigmahat-p(2))/p(1);

end
