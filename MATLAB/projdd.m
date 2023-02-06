function proj = projdd(params)
% Distance driven projection for equi-distance fanbeam CT
% Reference:
% B. De Man and S. Basu, "Distance-driven projection and
% backprojection in three dimensions," Physics in Medicine
% and Biology, vol. 49, pp. 2463-2475, 2004.
%
% input field of params:
%           'im'                         image
%           'Dsource2centr'     distance from source to center, mm
%           'Dsource2detec'    distance from source to detector, mm
%           'NumofBin'            number of detector bin
%           'pixelsize'              pixel size, mm
%           'binsize'                 detector bin size, mm
%           'binshift'                detector bin shift, mm
%           'iViews'                  view angels, degree
%    or   'NumofView'         number of projection view
%
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-10-09

%% fanbeam parameters

im = params.im;                                                       %im=phantom(512);
reconsize = size(im,1);                                               %reconsize=size(im,1); 
Dsource2centr = params.Dsource2centr;                                 %Dsource2centr=132;
Dsource2detec = params.Dsource2detec;                                 %Dsource2detec=180;
NumofBin = params.NumofBin;       % number of detector bins           %NumofBin=512;
pixelsize = params.pixelsize;                                         %pixelsize=0.075;
binsize = params.binsize;                                             %binsize=0.1;
binshift = params.binshift;               % detector shift, mm        %binshift=0;
                                                                      %NumofView=720;
% CT scanning views
if isfield(params, 'iViews')
    NumofView = length(params.iViews);
    iViews = params.iViews;
elseif isfield(params, 'NumofView')
    NumofView = params.NumofView;
    startAngle = 0;
    endAngle = 360-360/NumofView;
    iViews = startAngle : (endAngle - startAngle)/(NumofView - 1) : endAngle;               % degree
else
    error('There is no field ''iViews'' or ''NumofView'' ')
end

% if CT scanning is clockwise
if isfield(params, 'clockwise')
    if params.clockwise == 1
        iViews = mod(360 - iViews, 360);
    end
end

iBinPos = binsize*(-NumofBin/2:NumofBin/2) + binshift;    % rays' boudaries, NumofBin+1, mm
proj = zeros(NumofBin, NumofView);

%% main, equi-distance detector

for ip =  1:length(iViews)
    ip;
    iview = iViews(ip);
    ibeta = mod(iview/180*pi, 2*pi);
    
    % position of x-ray source
    xSource = -Dsource2centr*sin(ibeta);
    ySource = Dsource2centr*cos(ibeta);
    
    % projection for current view
    projline = zeros(NumofBin, 1);
    
    % source near y axis, so compute along x axis
    if (ibeta < pi/4) || (abs(ibeta - pi) < pi/4) || (2*pi - ibeta< pi/4)%makhroote bala o paeen
        
        % map all the rays' boundaries to x axis
        gammas = atan(iBinPos/Dsource2detec) ;
        detecPos = xSource + ySource*tan(gammas + ibeta);%%%Mori? mahalle barkhorde "khatte az source be marze detectorha" ba mehvare x
        if  abs(ibeta - pi) < pi/4%makhroote paeen
            % in this case, the position of detector is decreasing, we need
            % to reverse it
            detecPos = detecPos(end:-1:1);%%%Mori?
        end
        deteclen = detecPos(2:end) - detecPos(1:end-1);
        % cos of projection line
        cosRay = ySource./sqrt((detecPos - xSource).^2 + ySource^2);%%%Mori? Toole motefaavete ray ha dar tasvire morabaee
        
        % map all the image pixel boundaries to x axis
        linebound = ([0:reconsize] - reconsize/2);
        imbound = repmat(linebound, reconsize, 1);%mori x positions of boundaries of pixels
        imy = repmat([1:reconsize]' - reconsize/2 - 0.5, 1, reconsize+1); % y pos for each line pixels
        imPos = xSource + ySource*(imbound*pixelsize - xSource)./(ySource - imy*pixelsize);
        %mahalle barkhorde "khatte az source be noghate pixelha(sheklo negah konam)" ba mehvare x
        for ypix = 1:reconsize
            
            xPosLine = imPos(ypix,:) ;
            imline   =    im(ypix,:) ;
            
            currentPos = max([detecPos(1), xPosLine(1)]); % find the start position
            
            ix = 2;
            id = 2;           
            dPos = detecPos(id);
            xPos = xPosLine(ix);
            
            
            % if xPos <= currentPos, find start postion of image
            while xPos <= currentPos
                ix = ix + 1;
                if ix > reconsize+1
                    break;
                else
                    xPos = xPosLine(ix);
                end
            end
            
            % if dPos <= currentPos, find start position of detector
            while dPos <= currentPos
                id = id + 1;
                if id > NumofBin+1
                    break;
                else
                    dPos = detecPos(id);
                end
            end
            
            % main loop
            while ix <= reconsize + 1 && id <= NumofBin + 1
                xPos = xPosLine(ix);
                dPos = detecPos(id);
                id_1 = id - 1;
                ix_1 = ix - 1;
                if xPos < dPos
                    len = xPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(ix_1);
                    currentPos = xPos;
                    ix = ix + 1;
                elseif xPos > dPos
                    len = dPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(ix_1);
                    currentPos = dPos;
                    id = id + 1;
                else
                    len = xPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(ix_1);
                    currentPos = xPos;
                    ix = ix + 1;
                    id = id + 1;
                end
            end
            
        end
        
        
    else
        % map all the rays' boundaries to y axis
        gammas = atan(iBinPos/Dsource2detec) ;
        detecPos = ySource + xSource*tan(pi/2-(gammas + ibeta));
        if  abs(ibeta - 3/2*pi) <= pi/4
            % in this case, the position of detector is decreasing, we need
            % to reverse it
            detecPos = detecPos(end:-1:1);
        end
        deteclen = detecPos(2:end) - detecPos(1:end-1);
        % cos of projection line
        cosRay = xSource./sqrt((detecPos - ySource).^2 + xSource^2);
        
        % map all the image pixel boundaries to y axis
        linebound = ([0:reconsize] - reconsize/2)';
        imbound = repmat(linebound, 1, reconsize);
        imx = repmat([1:reconsize] - reconsize/2 - 0.5, reconsize+1, 1); % x pos for each line pixels
        imPos = ySource + xSource*(imbound*pixelsize - ySource)./(xSource - imx*pixelsize);
        
        for xpix = 1:reconsize
            
            iy = 2;
            id = 2;
            yPosLine = imPos(:, xpix);
            currentPos = max([detecPos(1), yPosLine(1)]); % find the start position
            dPos = detecPos(id);
            yPos = yPosLine(iy);
            imline = im(:, xpix);
            
            % if xPos <= currentPos, find start postion of image
            while yPos <= currentPos
                iy = iy + 1;
                if iy > reconsize+1
                    break;
                else
                    yPos = yPosLine(iy);
                end
            end
            
            % if dPos <= currentPos, find start position of detector
            while dPos <= currentPos
                id = id + 1;
                if id > NumofBin+1
                    continue;
                else
                    dPos = detecPos(id);
                end
            end
            
            % main loop
            while iy <= reconsize + 1 && id <= NumofBin + 1
                yPos = yPosLine(iy);
                dPos = detecPos(id);
                id_1 = id - 1;
                iy_1 = iy - 1;
                if yPos < dPos
                    len = yPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(iy_1);
                    currentPos = yPos;
                    iy = iy + 1;
                elseif yPos > dPos
                    len = dPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(iy_1);
                    currentPos = dPos;
                    id = id + 1;
                else
                    len = yPos - currentPos;
                    projline(id_1) =  projline(id_1) + len/deteclen(id_1)*imline(iy_1);
                    currentPos = yPos;
                    iy = iy + 1;
                    id = id + 1;
                end
            end
            
        end
        
    end
    
    if  abs(ibeta - pi) < pi/4 || abs(ibeta - 3/2*pi) <= pi/4
        % in this case, the position of detector is decreasing, we need
        % to reverse it
        projline = projline(end:-1:1);
        cosRay = cosRay(end:-1:1);
    end
    
    projline = pixelsize*projline./abs((cosRay(1:end-1)+cosRay(2:end))/2)';
    proj(:,ip) = projline;
    
end

% Tested OK
% def projdd(params):
%     pi = 3.14159265358979
%     #     % Distance driven projection for equi-distance fanbeam CT
%     #     % Reference:
%     #     % B. De Man and S. Basu, "Distance-driven projection and
%     #     % backprojection in three dimensions," Physics in Medicine
%     #     % and Biology, vol. 49, pp. 2463-2475, 2004.
%     #     %
%     #     % input field of params:
%     #     %           'im'                         image
%     #     %           'Dsource2centr'     distance from source to center, mm
%     #     %           'Dsource2detec'    distance from source to detector, mm
%     #     %           'NumofBin'            number of detector bin
%     #     %           'pixelsize'              pixel size, mm
%     #     %           'binsize'                 detector bin size, mm
%     #     %           'binshift'                detector bin shift, mm
%     #     %           'iViews'                  view angels, degree
%     #     %    or   'NumofView'         number of projection view
%     #     % Morteza Salehjahromi
%     im = params['im']
%     reconsize = im.shape[0]  
%     Dsource2centr = params['Dsource2centr']                            
%     Dsource2detec = params['Dsource2detec']                          
%     NumofBin = params['NumofBin']       
%     pixelsize = params['pixelsize']                                       
%     binsize = params['binsize']                                             
%     binshift = params['binshift']                                                                                  
%     # CT scanning views
%     if 'iViews' in params:
%         NumofView = params['iViews'].shape[0]
%         iViews = params['iViews']
%     elif 'NumofView' in params:
%         NumofView = params['NumofView']
%         startAngle = 0
%         endAngle = 360-360/NumofView;
%         iViews =  np.arange(startAngle,endAngle+1e-5,(endAngle - startAngle)/(NumofView - 1))
%     else:
%         raise ValueError('There is no field iViews or NumofView')
%     # if CT scanning is clockwise
%     if 'clockwise' in params:
%         if params['clockwise'] == 1:
%             iViews = np.mod(360 - iViews, 360)
%     # Not sure
%     iBinPos = binsize*np.arange(-NumofBin/2, NumofBin/2+1) + binshift    # rays' boudaries, NumofBin+1, mm
%     proj = np.zeros((NumofBin, NumofView))
%     ## main, equi-distance detector
%     for ip in range(0, iViews.size ):  #for ip =  1:length(iViews)
%         print(ip)
%         iview = iViews[ip]
%         ibeta = np.mod(iview/180*pi, 2*pi)
% 
%         # position of x-ray source
%         xSource = -Dsource2centr*np.sin(ibeta)
%         ySource = Dsource2centr*np.cos(ibeta)
% 
%         # projection for current view
%         projline = np.zeros((NumofBin,))
% 
%         # source near y axis, so compute along x axis
%         if (ibeta < pi/4) or (np.abs(ibeta - pi) < pi/4) or ((2*pi - ibeta)< pi/4): #makhroote bala o paeen
%             # map all the rays' boundaries to x axis
%             gammas = np.arctan(iBinPos/Dsource2detec)
%             detecPos = xSource + ySource*np.tan(gammas + ibeta);#Mori? mahalle barkhorde "khatte az source be marze detectorha" ba mehvare x
% 
%             if  np.abs(ibeta - pi) < pi/4:  #makhroote paeen
%                 # in this case, the position of detector is decreasing, we need
%                 # to reverse it
%                 detecPos = detecPos[::-1] #detecPos(end:-1:1)        #Mori?
%             deteclen = detecPos[1:] - detecPos[:-1] # detecPos(2:end) - detecPos(1:end-1);
%             # cos of projection line
%             cosRay = ySource/np.sqrt(np.power(detecPos - xSource,2) + np.power( ySource,2))   #Mori? Toole motefaavete ray ha dar tasvire morabaee
%             # map all the image pixel boundaries to x axis
%             linebound = np.arange(0, reconsize+1) - reconsize/2    # (1,513)  ([0:reconsize] - reconsize/2);
%             linebound = np.reshape(linebound,(1,linebound.shape[0]))
%             imbound = np.tile(linebound, (reconsize, 1)) #(512,513) #  repmat(linebound, reconsize, 1)     #mori x positions of boundaries of pixels
%             imy_temp = np.arange(1,reconsize+1) - reconsize/2 - 0.5
%             imy_temp = np.reshape(imy_temp,(imy_temp.shape[0],1))
%             imy = np.tile(imy_temp, (1, reconsize+1));
%             #imy = repmat([1:reconsize]' - reconsize/2 - 0.5, 1, reconsize+1); # y pos for each line pixels  (512,513)
%             imPos = xSource + ySource*(imbound*pixelsize - xSource)/(ySource - imy*pixelsize)
%             #mahalle barkhorde "khatte az source be noghate pixelha(sheklo negah konam)" ba mehvare x
%             for ypix in range(0,reconsize): 
%                 xPosLine = imPos[ypix,:] 
%                 imline   = im[ypix,:]
%                 currentPos = np.max([detecPos[0], xPosLine[0]]) # find the start position
%                 ix = 1;
%                 idd = 1;           
%                 dPos = detecPos[idd]
%                 xPos = xPosLine[ix]
%                 # if xPos <= currentPos, find start postion of image
%                 while xPos <= currentPos:
%                     ix = ix + 1;
%                     if ix > reconsize+1:
%                         break
%                     else:
%                         xPos = xPosLine[ix]
%                 # if dPos <= currentPos, find start position of detector
%                 while dPos <= currentPos:
%                     idd = idd + 1;
%                     if idd > NumofBin+1 :
%                         break
%                     else:
%                         dPos = detecPos[idd]
%                 # main loop
%                 while (ix <= reconsize ) and (idd <= NumofBin ):
%                     xPos = xPosLine[ix]
%                     dPos = detecPos[idd]
%                     id_1 = idd - 1
%                     ix_1 = ix - 1
%                     if xPos < dPos:
%                         lenn = xPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
%                         currentPos = xPos
%                         ix = ix + 1
%                     elif xPos > dPos:
%                         lenn = dPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
%                         currentPos = dPos
%                         idd = idd + 1
%                     else:
%                         lenn = xPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[ix_1]
%                         currentPos = xPos
%                         ix = ix + 1
%                         idd = idd + 1
%         else:
%             gammas = np.arctan(iBinPos/Dsource2detec)
%             detecPos = ySource + xSource*np.tan(pi/2-(gammas + ibeta))
%             if  np.abs(ibeta - 3/2*pi) <= pi/4:
%                 detecPos = detecPos[::-1] #(end:-1:1)
%             deteclen = detecPos[1:] - detecPos[:-1]
%             cosRay = xSource/np.sqrt( np.power(detecPos - ySource,2) + xSource**2)
%             linebound = (np.arange(0,reconsize+1) - reconsize/2)  #513*1  #([0:reconsize] - reconsize/2)';
%             linebound = np.reshape(linebound,(linebound.shape[0],1))
%             imbound = np.tile(linebound, (1,reconsize)) #(513,512)       # repmat(linebound, 1, reconsize);
%             imx_temp = np.arange(1,reconsize+1) - reconsize/2 - 0.5
%             imx_temp = np.reshape(imx_temp,(1, imx_temp.shape[0]))
%             imx = np.tile(imx_temp,(reconsize+1, 1))   #(513,512)
%             imPos = ySource + xSource*(imbound*pixelsize - ySource)/(xSource - imx*pixelsize)
%             for xpix in range(0,reconsize):
%                 iy = 1;
%                 idd = 1;
%                 yPosLine = imPos[:, xpix]
%                 currentPos = np.max([detecPos[0], yPosLine[0]]) # find the start position
%                 dPos = detecPos[idd]
%                 yPos = yPosLine[iy]
%                 imline = im[:, xpix]
% 
%                 # if xPos <= currentPos, find start postion of image
%                 while yPos <= currentPos:
%                     iy = iy + 1
%                     if iy > reconsize+1:
%                         break
%                     else:
%                         yPos = yPosLine[iy]
% 
%                 # if dPos <= currentPos, find start position of detector
%                 while dPos <= currentPos:
%                     idd = idd + 1
%                     if idd > NumofBin+1:
%                         continue
%                     else:
%                         dPos = detecPos[idd]
%                 # main loop
%                 while (iy <= reconsize ) and (idd <= NumofBin ):
%                     yPos = yPosLine[iy]
%                     dPos = detecPos[idd]
%                     id_1 = idd - 1
%                     iy_1 = iy - 1
%                     if yPos < dPos:
%                         lenn = yPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
%                         currentPos = yPos
%                         iy = iy + 1
%                     elif yPos > dPos:
%                         lenn = dPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
%                         currentPos = dPos
%                         idd = idd + 1
%                     else:
%                         lenn = yPos - currentPos
%                         projline[id_1] =  projline[id_1] + lenn/deteclen[id_1]*imline[iy_1]
%                         currentPos = yPos
%                         iy = iy + 1
%                         idd = idd + 1
%         if  (np.abs(ibeta - pi) < pi/4) or (np.abs(ibeta - 3/2*pi) <= pi/4):
%             projline = projline[::-1]    #(end:-1:1)
%             cosRay = cosRay[::-1]        #(end:-1:1)
%         temp_projline_div = np.abs((cosRay[0:-1]+cosRay[1:])/2)
%         projline = pixelsize*projline/temp_projline_div
%         proj[:,ip] = projline
%     return proj