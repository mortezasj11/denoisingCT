% % clc      %L13   0.075
% clear
% %%%load('ProjEnergy_noisyOrig_Simple5e3_8.mat')
% load('ProjEnergy_noisyOrig_Simple7e3_8.mat')
% NumofView=640;%NumofView=size(ProjEnergy,2);
% Proj1 = reshape(fliplr(rot90(ProjEnergy_noisy_Simple7e3_8(:,:,8),2)),512,1,NumofView);
% %%%%%%Proj1 = reshape(fliplr(rot90(ProjEnergy_Orig_Simple5e3_8(:,:,1),2)),512,1,NumofView);
% 
% Dsource2centr = 132;
% R=Dsource2centr;
% d = params.Dsource2detec-Dsource2centr;
% binsize   = 0.1;
% du=binsize;dv=binsize;
% pixelsize = 0.060 ;
% Reconsize = 512;
% [RecoResult1]=MortezaEquiDistanceSimulatedData(Proj1,R,d,du,dv,pixelsize,Reconsize);imshow(RecoResult1,[])
% 
% Ch8_8NoisySimple7e3=RecoResult1;
% save('Ch8_8NoisySimple7e3.mat','Ch8_8NoisySimple7e3')

function [RecoResult1]=MortezaEquiDistanceSimulatedData(Proj1,R,d,du,dv,pixelsize,Reconsize)
 %R=158;  d=255-158;  du=2*0.055; dv= 2*0.055; pixelsize = 18.41/512 ;  Reconsize=512;
% input:
%   dectectorshift: ȡֵ��Χ{0��1}����ⵥԪ���򣬼���ⵥԪ�����Ҹ����������������ڲ�����ģ��Forbildģ���ֵΪ1�����ڶ���ʵ�����ݸ�ֵΪ0

% % % % % %the later four para are [-1,1,-1,1] is the noraml mode,and means [-1,+1]*256 pixel and [-1,+1]*128 real size(mm);
YMin=-1;
YMax=+1;
XMin=-1;
XMax=+1;
shortScan=360;
Proj_Z=0;
%%%%%%%%%%%070827adding:pixelsize,for microCT reco,whose R is far smaller than the cxmax=512

[NumofRou, NumofRow, NumofView]=size(Proj1);
HalfNumofRou = round(NumofRou/2);%256
HalfNumofRow = round(NumofRow/2);%1
%HalfNumofView = round(NumofView/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters of Geometry Configuration
pi=3.14159265358979;
UnitofAngle=pi/180;

D=R+d;%D(255)     %R=DofSourceCentr(158)     d=Dobj2detec(97)
% % add short scan variables
% delta = atan(du*(NumofRou-1)/2/(d+R));                 % �ſ��Ƚǵ�һ��
% overlapbeta = shortScan*UnitofAngle - pi - 2*delta;    % ��ɨ���д�����С�걸ɨ��ĽǶ��� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ramp filter by FFT
B=1/2;%1/2*0.8; 
H=zeros(1,NumofRou*4);
for i=1:NumofRou*4
    n=i-NumofRou*2;
    %Shepp-Logan
%      H(i)=0.5*(4*B/pi)^2*(1-4*B*n*sin(pi/2*4*B*n))/(1-(4*B*n)^2);
    %Rama-Lak
    H(i)=2*B^2*sinc(2*n*B)-B^2*sinc(n*B)^2;
end
H=fshift(H,NumofRou*2+1);
H=fft(H);% dar vaghe in H, H*T^2 hast; choon T^2 emal nashode bood.
% figure
% plot(1:NumofRou*4,H);
% pause;

%%%%% Hilbert Transform  %%%%%
% du=dudown = dectsize * downsamp = 0.11 * 1 = 0.11
% u = u/alpha    ;  u   *  (1/alpha)= u * (158/256)
uline=((1:NumofRou) - HalfNumofRou)*du*(R/D); % du = [-255:256] * 0.11 * (158/255)
MhilbertProj1=zeros(2*NumofRou-1,NumofRow,NumofView);%(1023,1,360)
sinoProj1=zeros(1,NumofRou*4);%(1,2048)
position=(NumofRou/2)*3+1:(NumofRou/2)*5;   %769:1280  ;1280-769+1=512


for j=1:NumofRow
%     xi=(j-HalfNumofRow)*dv*(R/D);
    xi=Proj_Z*dv*(R/D);  %Proj_Z = 0 ; dv=0.11;  (%%%%%% ME: tooye jozve v tilt)
    Ruline=R./(R^2+uline.^2+xi^2).^0.5;  %size(uline)=512  
    for i=1:NumofView     %1:360
%         weight=ones(1,NumofRou);    % Parker weighting, �μ�S. Wesarg, M. Ebert, and T. Bortfeld, ��Parker weights revisited,��
%                                     % Medical Physics, vol. 29, no. 3, pp. 372-378, 2002
%         
%         alfa=[1:NumofRou]-HalfNumofRou;
%         alfa=(-1)^dectectorshift*atan(alfa*du/D);   
%         betaangle=(i-1)*shortScan/(NumofView-1)*UnitofAngle;
%         for m=1:NumofRou
%             
%             beta1 = betaangle/(2*delta - 2*alfa(m) + overlapbeta) - 0.5;
%             beta2 = (betaangle - pi + 2*alfa(m))/(2*delta + 2*alfa(m) + overlapbeta) - 0.5;
%             weight(m) = Sfunc(beta1) - Sfunc(beta2);
%             
%         end
        
        sinoProj1(position) = Ruline.*Proj1(:,j,i)'; %(1*512) .* (512*1)'
        sinofftProj1 = fft(sinoProj1);%1*2048
        
        MtempProj1 = real(ifft(sinofftProj1.*H)); %1*2048
        MhilbertProj1(:,j,i) = MtempProj1(NumofRou*1+2:NumofRou*3);%1023*1*360
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear H uline MtempProj1  sinoProj1;
% figure;
% imshow(reshape(MhilbertMuti,2*NumofRou-1,NumofView),[]);
% figure;
% imshow(reshape(MhilbertMono,2*NumofRou-1,NumofView),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBP of Generalized Feldkamp
cxmin=1;cxmax=Reconsize;%1,512, cxmax=reconsize
cymin=cxmin;cymax=cxmax;
czmin=180;czmax=180;
halfcx=round(cxmax/2);%256
OffX=round(cxmax/4*(XMax+XMin));%0
KX=(XMax-XMin)/2;%1
halfcy=halfcx;%256
OffY=round(cymax/4*(YMax+YMin));%0
KY=(YMax-YMin)/2;%1
IProj1=zeros(cxmax-cxmin+1,cymax-cymin+1,czmax-czmin+1);%512*512*1 or 512*512


% for beta=1:360*3
for beta=1:NumofView     %1.for each projection 
    beta
%     betaangle=(beta-1)/(NumofView/shortScan)*UnitofAngle;
%     shortScan=360;
    betaangle=(beta)/(NumofView/shortScan)*UnitofAngle;% 1:360 * (pi/180)

    zwave=0;    %zwave=z-h;    %zwave=120*p*UnitofAngle;    %zwave=-120*p*UnitofAngle;
     
    numz=0;
    for k=czmin:czmax
        numz=numz+1;

        numx=0;
        for i=cxmin:cxmax    %1:512
            numx=numx+1;
            x=((i-halfcx)*KX+OffX)*pixelsize;% x = (i-256)*pixelsize
            numy=0;
            for j=cymin:cymax   %1:512
                numy=numy+1;
                y=((j-halfcy)*KY+OffY)*pixelsize;  
                
                t=x*cos(betaangle)-y*sin(betaangle);
                s=x*sin(betaangle)+y*cos(betaangle);
                xi=R*zwave/(R-s);  %0
                v=round(xi*(D/R)/dv+HalfNumofRow); % dv=0.11, v=1;
                 if (v>0)&&(v<=NumofRow)
                    position=NumofRou+R*t/(R-s)*(D/R)/du;
                    positionup=ceil(position);
                    positionlow=floor(position);
                    if (positionup<2*NumofRou)&&(positionlow>0)
                        if (positionup~=position)
%                             InterResultProj1=MhilbertProj1(positionup,v,beta)*(position-positionlow)+MhilbertProj1(positionlow,v,beta)*(positionup-position);
                            InterResultProj1=MhilbertProj1(positionup,1,beta)*(position-positionlow)+MhilbertProj1(positionlow,1,beta)*(positionup-position);
                        else
%                             InterResultProj1=MhilbertProj1(position,v,beta);      
                            InterResultProj1=MhilbertProj1(position,1,beta);  
                        end
                            
%                     if position<2*NumofRou
%                         InterResult=Mhilbert(round(position),v,beta);               
                        IProj1(numx,numy,numz)=IProj1(numx,numy,numz)+InterResultProj1*(R/(R-s))^2;
                    end
                 end
                
            end%for k=1:cy           
        end%for j=1:cx
    end%for i=1:cz
end%for beta=1:360 
RecoResult1=IProj1/2*(2*pi/NumofView*shortScan/360)*(1/(du*R/D)); % ע�����ԭ����ȫɨ��IProj1Ҫ����2�����ڵ�Чֻ������pi�ڵ���Ч���ݲ��ó���2
RecoResult1=RecoResult1*10;     % ��λת����(1/cm)
% RecoResult1=upsidedown(RecoResult1);
%imshow(RecoResult1,[])
end
% horizontal vector is shifted towards right
function [M]=fshift(M,shiftnum)

    n=size(M);
    %row=n(1);
    col=n(2);

    if (shiftnum<1)||(shiftnum>=col)
        display('shiftnum is wrong!');
        return;
    end

    Mtemp=zeros(n);
    for i=1:col-shiftnum
        Mtemp(i+shiftnum)=M(i);
    end
    for i=(col-shiftnum+1):col
        Mtemp(i+shiftnum-col)=M(i);
    end

    M=Mtemp;
end
