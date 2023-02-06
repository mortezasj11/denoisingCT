import os
import  numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from scipy.io import loadmat
#from skimage.restoration import denoise_tv_chambolle,denoise_nl_means
def projdd(params):
    pi = 3.14159265358979
    #     % Distance driven projection for equi-distance fanbeam CT
    #     % Reference:
    #     % B. De Man and S. Basu, "Distance-driven projection and
    #     % backprojection in three dimensions," Physics in Medicine
    #     % and Biology, vol. 49, pp. 2463-2475, 2004.
    #     %
    #     % input field of params:
    #     %           'im'                         image
    #     %           'Dsource2centr'     distance from source to center, mm
    #     %           'Dsource2detec'    distance from source to detector, mm
    #     %           'NumofBin'            number of detector bin
    #     %           'pixelsize'              pixel size, mm
    #     %           'binsize'                 detector bin size, mm
    #     %           'binshift'                detector bin shift, mm
    #     %           'iViews'                  view angels, degree
    #     %    or   'NumofView'         number of projection view
    #     % Morteza Salehjahromi
    im = params['im']
    reconsize = im.shape[0]  
    Dsource2centr = params['Dsource2centr']                            
    Dsource2detec = params['Dsource2detec']                          
    NumofBin = params['NumofBin']       
    pixelsize = params['pixelsize']                                       
    binsize = params['binsize']                                             
    binshift = params['binshift']                                                                                  
    # CT scanning views
    if 'iViews' in params:
        NumofView = params['iViews'].shape[0]
        iViews = params['iViews']
    elif 'NumofView' in params:
        NumofView = params['NumofView']
        startAngle = 0
        endAngle = 360-360/NumofView;
        iViews =  np.arange(startAngle,endAngle+1e-5,(endAngle - startAngle)/(NumofView - 1))
    else:
        raise ValueError('There is no field iViews or NumofView')
    # if CT scanning is clockwise
    if 'clockwise' in params:
        if params['clockwise'] == 1:
            iViews = np.mod(360 - iViews, 360)
    # Not sure
    iBinPos = binsize*np.arange(-NumofBin/2, NumofBin/2+1) + binshift    # rays' boudaries, NumofBin+1, mm
    proj = np.zeros((NumofBin, NumofView))
    ## main, equi-distance detector
    for ip in range(0, iViews.size ):  #for ip =  1:length(iViews)
        print(ip)
        iview = iViews[ip]
        ibeta = np.mod(iview/180*pi, 2*pi)

        # position of x-ray source
        xSource = -Dsource2centr*np.sin(ibeta)
        ySource = Dsource2centr*np.cos(ibeta)

        # projection for current view
        projline = np.zeros((NumofBin,))

        # source near y axis, so compute along x axis
        if (ibeta < pi/4) or (np.abs(ibeta - pi) < pi/4) or ((2*pi - ibeta)< pi/4): #makhroote bala o paeen
            # map all the rays' boundaries to x axis
            gammas = np.arctan(iBinPos/Dsource2detec)
            detecPos = xSource + ySource*np.tan(gammas + ibeta);#Mori? mahalle barkhorde "khatte az source be marze detectorha" ba mehvare x

            if  np.abs(ibeta - pi) < pi/4:  #makhroote paeen
                # in this case, the position of detector is decreasing, we need
                # to reverse it
                detecPos = detecPos[::-1] #detecPos(end:-1:1)        #Mori?
            deteclen = detecPos[1:] - detecPos[:-1] # detecPos(2:end) - detecPos(1:end-1);
            # cos of projection line
            cosRay = ySource/np.sqrt(np.power(detecPos - xSource,2) + np.power( ySource,2))   #Mori? Toole motefaavete ray ha dar tasvire morabaee
            # map all the image pixel boundaries to x axis
            linebound = np.arange(0, reconsize+1) - reconsize/2    # (1,513)  ([0:reconsize] - reconsize/2);
            linebound = np.reshape(linebound,(1,linebound.shape[0]))
            imbound = np.tile(linebound, (reconsize, 1)) #(512,513) #  repmat(linebound, reconsize, 1)     #mori x positions of boundaries of pixels
            imy_temp = np.arange(1,reconsize+1) - reconsize/2 - 0.5
            imy_temp = np.reshape(imy_temp,(imy_temp.shape[0],1))
            imy = np.tile(imy_temp, (1, reconsize+1));
            #imy = repmat([1:reconsize]' - reconsize/2 - 0.5, 1, reconsize+1); # y pos for each line pixels  (512,513)
            imPos = xSource + ySource*(imbound*pixelsize - xSource)/(ySource - imy*pixelsize)
            #mahalle barkhorde "khatte az source be noghate pixelha(sheklo negah konam)" ba mehvare x
            for ypix in range(0,reconsize): 
                xPosLine = imPos[ypix,:] 
                imline   = im[ypix,:]
                currentPos = np.max([detecPos[0], xPosLine[0]]) # find the start position
                ix = 1;
                idd = 1;           
                dPos = detecPos[idd]
                xPos = xPosLine[ix]
                # if xPos <= currentPos, find start postion of image
                while xPos <= currentPos:
                    ix = ix + 1;
                    if ix > reconsize+1:
                        break
                    else:
                        xPos = xPosLine[ix]
                # if dPos <= currentPos, find start position of detector
                while dPos <= currentPos:
                    idd = idd + 1;
                    if idd > NumofBin+1 :
                        break
                    else:
                        dPos = detecPos[idd]
                # main loop
                while (ix <= reconsize ) and (idd <= NumofBin ):
                    xPos = xPosLine[ix]
                    dPos = detecPos[idd]
                    id_1 = idd - 1
                    ix_1 = ix - 1
                    if xPos < dPos:
                        lenn = xPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
                        currentPos = xPos
                        ix = ix + 1
                    elif xPos > dPos:
                        lenn = dPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
                        currentPos = dPos
                        idd = idd + 1
                    else:
                        lenn = xPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
                        currentPos = xPos
                        ix = ix + 1
                        idd = idd + 1
        else:
            gammas = np.arctan(iBinPos/Dsource2detec)
            detecPos = ySource + xSource*np.tan(pi/2-(gammas + ibeta))
            if  np.abs(ibeta - 3/2*pi) <= pi/4:
                detecPos = detecPos[::-1] #(end:-1:1)
            deteclen = detecPos[1:] - detecPos[:-1]
            cosRay = xSource/np.sqrt( np.power(detecPos - ySource,2) + xSource**2)
            linebound = (np.arange(0,reconsize+1) - reconsize/2)  #513*1  #([0:reconsize] - reconsize/2)';
            linebound = np.reshape(linebound,(linebound.shape[0],1))
            imbound = np.tile(linebound, (1,reconsize)) #(513,512)       # repmat(linebound, 1, reconsize);
            imx_temp = np.arange(1,reconsize+1) - reconsize/2 - 0.5
            imx_temp = np.reshape(imx_temp,(1, imx_temp.shape[0]))
            imx = np.tile(imx_temp,(reconsize+1, 1))   #(513,512)
            imPos = ySource + xSource*(imbound*pixelsize - ySource)/(xSource - imx*pixelsize)
            for xpix in range(0,reconsize):
                iy = 1;
                idd = 1;
                yPosLine = imPos[:, xpix]
                currentPos = np.max([detecPos[0], yPosLine[0]]) # find the start position
                dPos = detecPos[idd]
                yPos = yPosLine[iy]
                imline = im[:, xpix]

                # if xPos <= currentPos, find start postion of image
                while yPos <= currentPos:
                    iy = iy + 1
                    if iy > reconsize+1:
                        break
                    else:
                        yPos = yPosLine[iy]

                # if dPos <= currentPos, find start position of detector
                while dPos <= currentPos:
                    idd = idd + 1
                    if idd > NumofBin+1:
                        continue
                    else:
                        dPos = detecPos[idd]
                # main loop
                while (iy <= reconsize ) and (idd <= NumofBin ):
                    yPos = yPosLine[iy]
                    dPos = detecPos[idd]
                    id_1 = idd - 1
                    iy_1 = iy - 1
                    if yPos < dPos:
                        lenn = yPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
                        currentPos = yPos
                        iy = iy + 1
                    elif yPos > dPos:
                        lenn = dPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
                        currentPos = dPos
                        idd = idd + 1
                    else:
                        lenn = yPos - currentPos
                        projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
                        currentPos = yPos
                        iy = iy + 1
                        idd = idd + 1
        if  (np.abs(ibeta - pi) < pi/4) or (np.abs(ibeta - 3/2*pi) <= pi/4):
            projline = projline[::-1]    #(end:-1:1)
            cosRay = cosRay[::-1]        #(end:-1:1)
        temp_projline_div = np.abs((cosRay[0:-1]+cosRay[1:])/2)
        projline = pixelsize*projline/temp_projline_div
        proj[:,ip] = projline
    return proj

def fshift(M,shiftnum):
    # horizontal vector is shifted towards right
    n = M.shape
    #row = n[0]
    col = n[1]
    if (shiftnum<1) or (shiftnum>=col):
        raise ValueError('shiftnum is wrong!')
        return
    Mtemp= np.zeros(n)
    for i in range(0,col-shiftnum):
        Mtemp[0,i+shiftnum]=M[0,i]
    for i in range(col-shiftnum,col):
        Mtemp[0,i+shiftnum-col]=M[0,i]
    return Mtemp

def MortezaEquiDistanceRealDATA(Proj1,R,d,du,dv,pixelsize,Reconsize):
    YMin = -1
    YMax = +1
    XMin = -1
    XMax = +1
    shortScan = 360
    Proj_Z = 0
    NumofRou, NumofRow, NumofView = Proj1.shape
    HalfNumofRou = np.round(NumofRou/2) #256
    HalfNumofRow = 1 # np.round(NumofRow/2) #1  ?????????????????????????????????????????
    pi = 3.14159265358979
    UnitofAngle = pi/180
    D = R+d     #D(255)     %R=DofSourceCentr(158)     d=Dobj2detec(97)
    # Ramp filter by FFT
    B = 1/2                 #1/2*0.8; 
    
    H = np.zeros((1,NumofRou*4))
    
    for i in range(0, NumofRou*4):  #=1:NumofRou*4
        n = i+1-NumofRou*2
        H[0,i] = 2*B**2*np.sinc(2*n*B)-B**2*np.sinc(n*B)**2
      
    H = fshift(H,NumofRou*2+1)
    
    H = np.fft.fft(H) # dar vaghe in H, H*T^2 hast; choon T^2 emal nashode bood.
    uline = (np.arange(1,NumofRou+1) - HalfNumofRou)*du*(R/D) # du = [-255:256] * 0.11 * (158/255)
    MhilbertProj1 = np.zeros((2*NumofRou-1,NumofRow,NumofView)) #(1023,1,360)
    sinoProj1 = np.zeros((1,NumofRou*4))#(1,2048)
    position = np.int16(np.arange( (NumofRou/2)*3+1 , (NumofRou/2)*5+1e-4 ))   #769:1280  ;1280-769+1=512
    
    for j in range(0,NumofRow):
        #xi=(j-HalfNumofRow)*dv*(R/D);
        xi = Proj_Z*dv*(R/D) #Proj_Z = 0 ; dv=0.11;  (%%%%%% ME: tooye jozve v tilt)
        Ruline=R/(R**2 + np.power(uline,2) + xi**2)**0.5  #size(uline)=512
        
        for i in range(0,NumofView):

            sinoProj1[0, position] = Ruline*Proj1[:,j,i]  #(1*512) .* (512*1)'
            sinofftProj1 = np.fft.fft(sinoProj1)         #1*2048

            MtempProj1 = np.real(np.fft.ifft(sinofftProj1*H))            #1*2048
            
            MhilbertProj1[:,j,i] = MtempProj1[0,NumofRou*1+2:NumofRou*3+1]   #1023*1*360   

    del H; del uline;del MtempProj1;del  sinoProj1;
    # FBP of Generalized Feldkamp
    cxmin = 1; cxmax= int(Reconsize)  #1,512, cxmax=reconsize
    cymin = cxmin ; cymax = cxmax
    czmin = 180 ; czmax = 180
    halfcx = np.round(cxmax/2) #256
    OffX = np.round(cxmax/4*(XMax+XMin)) #0
    KX = (XMax-XMin)/2 #1
    halfcy = halfcx #256
    OffY = np.round(cymax/4*(YMax+YMin))#0
    KY = (YMax-YMin)/2 #1
    IProj1 = np.zeros((int(cxmax-cxmin+1),int(cymax-cymin+1),int(czmax-czmin+1))) #512*512*1 or 512*512
    
    # for beta=1:360*3
    for beta in range(0, NumofView):    #1.for each projection 
        print(beta)
        betaangle=(beta+1)/(NumofView/shortScan)*UnitofAngle #1:360 * (pi/180)
        zwave = 0   #zwave=z-h;    %zwave=120*p*UnitofAngle;    %zwave=-120*p*UnitofAngle;
        numz  = -1
        
        for k in range (czmin,czmax+1):
            numz = numz+1
            numx = -1
            
            for i in range(cxmin,cxmax+1):  #1:512
                numx = numx + 1
                x = ((i-halfcx)*KX + OffX)*pixelsize # x = (i-256)*pixelsize
                numy = -1
                
                for j in range(cymin, cymax):   #1:512
                    numy=numy+1;
                    y=((j-halfcy)*KY+OffY)*pixelsize
                    t = x*np.cos(betaangle) - y*np.sin(betaangle)
                    s = x*np.sin(betaangle) + y*np.cos(betaangle)
                    xi = R*zwave/(R-s)  #0
                    v = np.round(xi*(D/R)/dv + HalfNumofRow) # dv=0.11, v=1;
                    #print(numy,t,s,xi,v)
                    #breakpoint()
                    if (v>0) and (v<=NumofRow):
                        position = NumofRou + R*t/(R-s)*(D/R)/du
                        positionup = int(np.ceil(position))
                        positionlow = int(np.floor(position))
                        #breakpoint()
                        if (positionup<2*NumofRou) and (positionlow>0):
                            if positionup != position:
                                InterResultProj1 = MhilbertProj1[positionup,0,beta]*(position-positionlow) + MhilbertProj1[positionlow,0,beta]*(positionup-position)
                            else:
                                
                                InterResultProj1 = MhilbertProj1[int(position),0,beta]       
                            
                            IProj1[numx,numy,numz] = IProj1[numx,numy,numz] + InterResultProj1*((R/(R-s))**2)

    RecoResult1 = IProj1/2*(2*pi/NumofView*shortScan/360)*(1/(du*R/D))
    RecoResult1 = RecoResult1    #(1/cm)
    return RecoResult1

# PARAMETERS 
params = {}
params['Dsource2centr'] = 165
params['Dsource2detec'] = 180
params['NumofView']   = 256
params['NumofBin']    = 800     # number of detector bins
params['pixelsize' ]  = 0.05
params['binsize' ]    = 0.05
params['binshift' ]   = 0       # detector shift, mm
coef = 1
params['pixelsize' ]  = 0.05*coef
params['reconsize']   = 512/coef


path = 'C:/Users/MSalehjahromi/ZipingRongqin/2_Matching_cohort/CT_machingOld'
ct4  = '868417_150814.nii.gz'
path_cti = os.path.join(path, ct4)
ct_w = nib.load(path_cti).get_fdata()
slice_n = ct_w.shape[2]//2
params['im'] = ct_w[:, :, slice_n]  # np.ones((512,512)) # tv2_nl # np.ones((512,512))
Proj = projdd(params);
plt.imshow(Proj) #plt.imshow(Proj,vmin=0,vmax=3.5)


#Proj = loadmat('Proj16.mat')['Proj']


##################################################################################
Proj1 = np.reshape(np.flipud(Proj),(params['NumofBin'],1,params['NumofView']) )
proj = np.rot90(Proj1,-1)

NumofRoudown = params['NumofBin']
Proj1 = np.reshape(np.fliplr(np.rot90(proj)),(NumofRoudown,1,params['NumofView']) )

R = params['Dsource2centr']
d = params['Dsource2detec'] - params['Dsource2centr']
du = params['binsize']
dv = params['binsize']

coef = 1
pixelsize = coef*params['pixelsize' ] 
Reconsize =512/coef
RecoResult1 = MortezaEquiDistanceRealDATA(Proj1,R,d,du,dv,pixelsize,Reconsize)

plt.imshow(np.rot90(RecoResult1,-1),vmin=-1000,vmax=300, cmap='gray')