
%                                Proj1 = reshape(flipud(ProjEnergy_noisyPix06_2e3(:,:,1)),512,1,NumofView);
%%Input proj of real data
%         load('ProjReal1_8.mat')
%         proj=rot90(ProjReal1_8(:,:,1),-1);%imshow(proj,[])
% NumofView=360;
% NumofRoudown=512;
% Proj1 = reshape(fliplr(rot90(proj)),NumofRoudown,1,NumofView);
% R=158;
% d=255-158;
% du=2*0.055;
% dv= 2*0.055;clc
% pixelsize = 18.41/512 ;
% Reconsize=512;
%[RecoResult1]=MortezaEquiDistanceRealDATA(Proj1,R,d,du,dv,pixelsize,Reconsize);imshow(rot90(RecoResult1,-1),[])
function [RecoResult1]=MortezaEquiDistanceRealDATA(Proj1,R,d,du,dv,pixelsize,Reconsize)
 %R=158;  d=255-158;  du=2*0.055; dv= 2*0.055; pixelsize = 18.41/512 ;  Reconsize=512;
% input:
%   dectectorshift: 取值范围{0，1}，检测单元反向，即检测单元左正右负还是左负右正，对于产生的模拟Forbild模体该值为1，对于多数实际数据该值为0

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
% delta = atan(du*(NumofRou-1)/2/(d+R));                 % 张开扇角的一半
% overlapbeta = shortScan*UnitofAngle - pi - 2*delta;    % 短扫描中大于最小完备扫描的角度数 
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
%         weight=ones(1,NumofRou);    % Parker weighting, 参见S. Wesarg, M. Ebert, and T. Bortfeld, “Parker weights revisited,”
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
%     betaangle=(beta-1)/(NumofView/shortScan)*UnitofAngle;
%     shortScan=360;
beta
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
RecoResult1=IProj1/2*(2*pi/NumofView*shortScan/360)*(1/(du*R/D)); % 注意对于原来的全扫描IProj1要除以2，现在等效只利用了pi内的有效数据不用除以2
RecoResult1=RecoResult1;     % 单位转换成(1/cm)
% RecoResult1=upsidedown(RecoResult1);
%imshow(RecoResult1,[])
CH_1_2e3=RecoResult1;
save('CH_1_2e3.mat','CH_1_2e3')




%# Tested OK

% def fshift(M,shiftnum):
%     # horizontal vector is shifted towards right
%     n = M.shape
%     #row = n[0]
%     col = n[1]
%     if (shiftnum<1) or (shiftnum>=col):
%         raise ValueError('shiftnum is wrong!')
%         return
%     Mtemp= np.zeros(n)
%     for i in range(0,col-shiftnum):
%         Mtemp[0,i+shiftnum]=M[0,i]
%     for i in range(col-shiftnum,col):
%         Mtemp[0,i+shiftnum-col]=M[0,i]
%     return Mtemp
% 
% def MortezaEquiDistanceRealDATA(Proj1,R,d,du,dv,pixelsize,Reconsize):
%     YMin = -1
%     YMax = +1
%     XMin = -1
%     XMax = +1
%     shortScan = 360
%     Proj_Z = 0
%     NumofRou, NumofRow, NumofView = Proj1.shape
%     HalfNumofRou = np.round(NumofRou/2) #256
%     HalfNumofRow = 1 # np.round(NumofRow/2) #1  ?????????????????????????????????????????
%     pi = 3.14159265358979
%     UnitofAngle = pi/180
%     D = R+d     #D(255)     %R=DofSourceCentr(158)     d=Dobj2detec(97)
%     # Ramp filter by FFT
%     B = 1/2                 #1/2*0.8; 
%     
%     H = np.zeros((1,NumofRou*4))
%     
%     for i in range(0, NumofRou*4):  #=1:NumofRou*4
%         n = i+1-NumofRou*2
%         H[0,i] = 2*B**2*np.sinc(2*n*B)-B**2*np.sinc(n*B)**2
%       
%     H = fshift(H,NumofRou*2+1)
%     
%     H = np.fft.fft(H) # dar vaghe in H, H*T^2 hast; choon T^2 emal nashode bood.
%     uline = (np.arange(1,NumofRou+1) - HalfNumofRou)*du*(R/D) # du = [-255:256] * 0.11 * (158/255)
%     MhilbertProj1 = np.zeros((2*NumofRou-1,NumofRow,NumofView)) #(1023,1,360)
%     sinoProj1 = np.zeros((1,NumofRou*4))#(1,2048)
%     position = np.int16(np.arange( (NumofRou/2)*3+1 , (NumofRou/2)*5+1e-4 ))   #769:1280  ;1280-769+1=512
%     
%     for j in range(0,NumofRow):
%         #xi=(j-HalfNumofRow)*dv*(R/D);
%         xi = Proj_Z*dv*(R/D) #Proj_Z = 0 ; dv=0.11;  (%%%%%% ME: tooye jozve v tilt)
%         Ruline=R/(R**2 + np.power(uline,2) + xi**2)**0.5  #size(uline)=512
%         
%         for i in range(0,NumofView):
% 
%             sinoProj1[0, position] = Ruline*Proj1[:,j,i]  #(1*512) .* (512*1)'
%             sinofftProj1 = np.fft.fft(sinoProj1)         #1*2048
% 
%             MtempProj1 = np.real(np.fft.ifft(sinofftProj1*H))            #1*2048
%             
%             MhilbertProj1[:,j,i] = MtempProj1[0,NumofRou*1+2:NumofRou*3+1]   #1023*1*360   
% 
%     del H; del uline;del MtempProj1;del  sinoProj1;
%     # FBP of Generalized Feldkamp
%     cxmin = 1; cxmax= int(Reconsize)  #1,512, cxmax=reconsize
%     cymin = cxmin ; cymax = cxmax
%     czmin = 180 ; czmax = 180
%     halfcx = np.round(cxmax/2) #256
%     OffX = np.round(cxmax/4*(XMax+XMin)) #0
%     KX = (XMax-XMin)/2 #1
%     halfcy = halfcx #256
%     OffY = np.round(cymax/4*(YMax+YMin))#0
%     KY = (YMax-YMin)/2 #1
%     IProj1 = np.zeros((int(cxmax-cxmin+1),int(cymax-cymin+1),int(czmax-czmin+1))) #512*512*1 or 512*512
%     
%     # for beta=1:360*3
%     for beta in range(0, NumofView):    #1.for each projection 
%         print(beta)
%         betaangle=(beta+1)/(NumofView/shortScan)*UnitofAngle #1:360 * (pi/180)
%         zwave = 0   #zwave=z-h;    %zwave=120*p*UnitofAngle;    %zwave=-120*p*UnitofAngle;
%         numz  = -1
%         
%         for k in range (czmin,czmax+1):
%             numz = numz+1
%             numx = -1
%             
%             for i in range(cxmin,cxmax+1):  #1:512
%                 numx = numx + 1
%                 x = ((i-halfcx)*KX + OffX)*pixelsize # x = (i-256)*pixelsize
%                 numy = -1
%                 
%                 for j in range(cymin, cymax):   #1:512
%                     numy=numy+1;
%                     y=((j-halfcy)*KY+OffY)*pixelsize
%                     t = x*np.cos(betaangle) - y*np.sin(betaangle)
%                     s = x*np.sin(betaangle) + y*np.cos(betaangle)
%                     xi = R*zwave/(R-s)  #0
%                     v = np.round(xi*(D/R)/dv + HalfNumofRow) # dv=0.11, v=1;
%                     #print(numy,t,s,xi,v)
%                     #breakpoint()
%                     if (v>0) and (v<=NumofRow):
%                         position = NumofRou + R*t/(R-s)*(D/R)/du
%                         positionup = int(np.ceil(position))
%                         positionlow = int(np.floor(position))
%                         #breakpoint()
%                         if (positionup<2*NumofRou) and (positionlow>0):
%                             if positionup != position:
%                                 InterResultProj1 = MhilbertProj1[positionup,0,beta]*(position-positionlow) + MhilbertProj1[positionlow,0,beta]*(positionup-position)
%                             else:
%                                 
%                                 InterResultProj1 = MhilbertProj1[int(position),0,beta]       
%                             
%                             IProj1[numx,numy,numz] = IProj1[numx,numy,numz] + InterResultProj1*((R/(R-s))**2)
% 
%     RecoResult1 = IProj1/2*(2*pi/NumofView*shortScan/360)*(1/(du*R/D))
%     RecoResult1 = RecoResult1    #(1/cm)
%     return RecoResult1
