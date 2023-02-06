% TNN1����tensor����ת����mode-3�����svd�ֽ⣬������ֵ���д���
%
% 2014-1-7

function imgkMCTNN1 = TNN1(imgkMC, nNormThresh, typeflag)

% typeflag�� ��־λ��ȡֵ'nNorm'��'Tresh'

% % ��ֵ��ʾԤ����������ֵ����֮�������ֵ���㣬ȡֵ��Χ��int [1, nMC-1]��
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

% �Էֽ�õ�������ֵ��������ֵ����
if strcmp(typeflag, 'nNorm')
    % ָ������������ֵ
    tol = S(nNormThresh+1, nNormThresh+1);
elseif strcmp(typeflag, 'Thresh')
    % ָ����ֵ
    tol = nNormThresh;
end
Ssoft = S - tol;
Ssoft(Ssoft < 0) = 0;
imgkMatrixNuclear = U*Ssoft*V';

% �������ľ���ָ���tensor������ʾ
imgkMCTNN1 = zeros(rows, cols, nMC);
for iMC = 1:nMC
    imgkMCTNN1(:,:,iMC) = reshape(imgkMatrixNuclear(iMC, :), rows, cols);
end


