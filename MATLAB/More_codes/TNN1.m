% TNN1，经tensor数据转化成mode-3矩阵后svd分解，对特征值进行处理
%
% 2014-1-7

function imgkMCTNN1 = TNN1(imgkMC, nNormThresh, typeflag)

% typeflag： 标志位，取值'nNorm'和'Tresh'

% % 该值表示预留几个特征值，这之后的特征值置零，取值范围：int [1, nMC-1]。
% nNorm = 1;

[rows, cols, nMC] = size(imgkMC);

if strcmp(typeflag, 'nNorm') && (nNormThresh <1 || nNormThresh>=nMC)
    disp('Error: Please input right nNormThresh !!!')
    return
end

imgkMatrix = zeros(nMC, rows*cols);
for iMC = 1:nMC
    imgkMatrix(iMC, :) = reshape(imgkMC(:,:,iMC), 1, rows*cols);
end
[U, S, V] = mySVD(imgkMatrix);
S = full(S);

% 对分解得到的奇异值进行软阈值处理
if strcmp(typeflag, 'nNorm')
    % 指定保留的奇异值
    tol = S(nNormThresh+1, nNormThresh+1);
elseif strcmp(typeflag, 'Thresh')
    % 指定阈值
    tol = nNormThresh;
end
Ssoft = S - tol;
Ssoft(Ssoft < 0) = 0;
imgkMatrixNuclear = U*Ssoft*V';

% 将处理后的矩阵恢复成tensor，并显示
imgkMCTNN1 = zeros(rows, cols, nMC);
for iMC = 1:nMC
    imgkMCTNN1(:,:,iMC) = reshape(imgkMatrixNuclear(iMC, :), rows, cols);
end


