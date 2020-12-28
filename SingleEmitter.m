%% Simulate Emitters
% v2.1 | based on SNR = mean(PL)/std(BG)
% Introduce more noise/variance between frames
% September 8, 2020
% by Wesley Chiang
% Target SNR QD=5.33 | R6G = 2.02
%% Parameter setup
ROI = 8;        % 8x8 pixel ROI equivalent to that used to get traces
ROI_PL = repmat(false,ROI,ROI); 
ROI_PL(2:end-1,2:end-1)=true;  % 6x6 region used to get PL traces
ROI_BG = ~ROI_PL;           % region wrapping around PL ROI to get BG trace
A = [sum(sum(ROI_PL)) sum(sum(ROI_BG))];    % A represents are of ROI
[xPL,yPL] = find(ROI_PL);   % coordinates for PL ROI
[xBG,yBG] = find(ROI_BG);   % coordinates for BG ROI
bits = 16;      % bit depth used for image
QE = 0.8;        % Quantum effiency of pco edge 4.2 @ 600nm = 80%
depth = 45000;  % Quantum Well Depth of pco edge 4.2 = 45000 e-
noise = 2.3;    % rms readout noise of pco edge 4.2 = 2.3 e-
frames = 1000;  % number of frames used in traces
fps = 40;       % frame rate of traces
%% Save location
F_save = uigetdir(pwd,'Set Save Location');
%% Create emitter of SNR and SNR 16

m = [1:ROI]-(ROI+1)/2;  % grid layout of length ROI centered at 0
[x,y] = meshgrid(m,m);  % m-by-m grid of values
ind = x.^2+y.^2;        % square of indices of size ROI
emit = ind<((ROI/4)^2);
ap = ind<((ROI/8)^2);
n = 1; figure(n); n=n+1; 
subplot(2,1,1); imagesc(emit); colormap(gray);axis image; title('Emitter');
subplot(2,1,2); imagesc(ap); colormap(gray); axis image; title('Aperture');

clim = [0 15000];       % Scale upper and lower bounds for plot

emitLow = 500*double(emit); % emitter low SNR
emitLow(~emit)=1; %emit8 = 500*emit8;

emitHigh = 350*double(emit);% emitter high SNR
emitHigh(~emit)=1; %emit3 = 500*emit3;
figure(n); n=n+1; 
subplot(2,1,1); imagesc(emitLow); title('Emitter Low SNR');
colormap(gray); colorbar; axis image
subplot(2,1,2); imagesc(emitHigh); title('Emitter High SNR');
colormap(gray); colorbar; axis image

%% Add blur to emitters
radcoeff = 2*pi/((ROI/16)^2); q_phase = exp(1i*radcoeff*ind);
    % Generate quadrat
h = fftshift(fftn(ifftshift(ap.*q_phase)));
PSF = abs(h).^2;
H = fftshift(fftn(ifftshift(PSF)));

F_Low = fftshift(fftn(ifftshift(emitLow)));
B_Low = F_Low.*H;
%b_Low = 5*fftshift(ifftn(ifftshift(B_Low)));
b_Low = fftshift(ifftn(ifftshift(B_Low)));
rLow = mean(b_Low(ROI_PL))/std(b_Low(ROI_BG));

F_High = fftshift(fftn(ifftshift(emitHigh)));
B_High = F_High.*H;
%b_High = 5*fftshift(ifftn(ifftshift(B_High)));
b_High = fftshift(ifftn(ifftshift(B_High)));
rHigh = mean(b_High(ROI_PL))/std(b_High(ROI_BG));

figure(n); n=n+1;
subplot(2,1,1); imagesc(1/4*b_Low);
rLow = mean(b_Low(ROI_PL))/std(b_Low(ROI_BG));
tit = sprintf('Blurred emitter SNR = %.2f',rLow);
title(tit);
colormap(gray); colorbar; axis image
subplot(2,1,2); imagesc(1/4*b_High); 
rHigh = mean(b_High(ROI_PL))/std(b_High(ROI_BG));
tit = sprintf('Blurred emitter SNR = %.2f', rHigh);
title(tit);
colormap(gray); colorbar; axis image
%% Create n-frame movies of emitters
movLow = repmat(1/4*b_Low,1,1,frames);
movHigh = repmat(1/2*b_High,1,1,frames);
%% Introduce desired SNR
% Low SNR
scaleBG = 20/rLow*(-0.5+(1.25+0.5)*rand(size(movLow)));
movLow = movLow.*ROI_BG.*scaleBG+movLow; movLow(movLow<0) = 0;
    % Increase noise in BG to achieve desired SNR
PL_Low = mean(reshape(movLow(repmat(ROI_PL,1,1,frames)),A(1),frames));
BG_Low = mean(reshape(movLow(repmat(ROI_BG,1,1,frames)),A(2),frames));
BGSD_Low = std(reshape(movLow(repmat(ROI_BG,1,1,frames)),A(2),frames));
s2n_Low = mean(PL_Low./BGSD_Low);
scalePL = sqrt(s2n_Low)*(0.35+(1-0.35)*rand(size(movLow)));
movLow = 2/3*movLow.*ROI_PL.*scalePL+movLow.*ROI_BG;
    % Introduce random noise across emitter pixels to create frame to frame
    % variance that ~SNR
PL_Low = mean(reshape(movLow(repmat(ROI_PL,1,1,frames)),A(1),frames));
s2n_Low = mean(PL_Low./BGSD_Low);

% High SNR
scaleBG = 4/rHigh*(-1.5+(1.75+1.5)*rand(size(movHigh)));
movHigh = movHigh.*ROI_BG.*scaleBG+movHigh;
    % Noise in BG
PL_High = mean(reshape(movHigh(repmat(ROI_PL,1,1,frames)),A(1),frames));
BG_High = mean(reshape(movHigh(repmat(ROI_BG,1,1,frames)),A(2),frames));
BGSD_High = std(reshape(movHigh(repmat(ROI_BG,1,1,frames)),A(2),frames));
s2n_High = mean(PL_High./BGSD_High);
scalePL = 1/2*sqrt(s2n_High)*(0.4+(1-0.4)*rand(size(movHigh)));
movHigh = movHigh.*ROI_PL.*scalePL+movHigh.*ROI_BG;
    % Frame to frame variance
PL_High = mean(reshape(movHigh(repmat(ROI_PL,1,1,frames)),A(1),frames));
s2n_High = mean(PL_High./BGSD_High);

% Plot
figure(n); n=n+1;
subplot(2,1,1); plot(PL_Low,'r'); hold on; plot(BG_Low,'k');
tit = sprintf('No Blink Trace for SNR = %.2f',s2n_Low);
title(tit);
subplot(2,1,2); plot(PL_High,'r'); hold on; plot(BG_High,'k');
tit = sprintf('No Blink Trace for SNR = %.2f',s2n_High);
title(tit);

%% Introduce Poisson-Gaussian Noise due to camera parameters
PLow = poissrnd(round(movLow*QE*depth/(2^bits)));
PHigh = poissrnd(round(movHigh*QE*depth/(2^bits)));
    % Poisson (shot noise)
    % Random poisson noising ontop of photoelectrons from simulated emitter
    % Quantum efficiency and well depth of pco edge 4.2
GPLow = PLow+sqrt(noise)*randn(size(PLow)); 
GPHigh = PHigh+sqrt(noise)*randn(size(PHigh));
    % Gaussian + Poisson noise
    
figure(n); n=n+1;
subplot(2,3,1); imagesc(GPLow(:,:,1)); colorbar; axis image;  
title('GP Low SNR - Frame 1');
subplot(2,3,2); imagesc(GPLow(:,:,300)); colorbar; axis image;  
title('GP Low SNR - Frame 300');
subplot(2,3,3); imagesc(GPLow(:,:,800)); colorbar; axis image;  
title('GP Low SNR - Frame 800');
subplot(2,3,4); imagesc(GPHigh(:,:,1)); colorbar; axis image;  
title('GP High SNR - Frame 1');
subplot(2,3,5); imagesc(GPHigh(:,:,300)); colorbar; axis image;  
title('GP High SNR  - Frame 300');
subplot(2,3,6); imagesc(GPHigh(:,:,800)); colorbar; axis image;  
title('GP High SNR  - Frame 800');
colormap(gray);

%% Blink Parameters
bps = [2,2.7];                  % blinks per second
numOff = round(bps/40*frames);  % number of frames to be off

%% 2 blinks per second
offFrames = round(frames*rand(1,numOff(1))); % frames off @ 2bps
blinkLow = GPLow; blinkHigh = GPHigh;
BG_Low = mean(reshape(blinkLow(repmat(ROI_BG,1,1,frames)),A(2),frames));
BG_High = mean(reshape(blinkHigh(repmat(ROI_BG,1,1,frames)),A(2),frames));
BGSD_Low = std(reshape(blinkLow(repmat(ROI_BG,1,1,frames)),A(2),frames));
BGSD_High = std(reshape(blinkHigh(repmat(ROI_BG,1,1,frames)),A(2),frames));
    % reshape = A-by-frames array where each column represents a frame
    % mean = 1-by-frames array where each column is the grey level of frame
offLow = mean(BG_Low)+sqrt(mean(BGSD_Low))*randn(size(ROI_PL)); 
offHigh = mean(BG_High)+sqrt(mean(BGSD_High))*randn(size(ROI_PL));
offLow(offLow<0)=0; offHigh(offHigh<0)=0;
blinkLow(:,:,offFrames) = repmat(offLow,1,1,numOff(1));
    % 2bps @ SNR L
blinkHigh(:,:,offFrames) = repmat(offHigh,1,1,numOff(1));
    % 2bps @ SNR H
    % Set coordinates of PL ROI to the average BG level @ SNR


fLow = fullfile(F_save,'LowSNR_2bps.ome.tiff');
bfsave(uint16(blinkLow),fLow);
fHigh = fullfile(F_save,'HighSNR_2bps.ome.tiff');
bfsave(uint16(blinkHigh),fHigh);
    %}

%% Check 2bps trace
PL_Low = mean(reshape(blinkLow(repmat(ROI_PL,1,1,frames)),A(1),frames));
PL_High = mean(reshape(blinkHigh(repmat(ROI_PL,1,1,frames)),A(1),frames));
    % reshape = A-by-frames array where each column represents a frame
    % mean = 1-by-frames array where each column is the grey level of frame    
s2n = [mean(PL_Low./BGSD_Low) mean(PL_High./BGSD_High)]; % SNR of the traces
figure(n); n=n+1;
subplot(2,1,1); tit = sprintf('Trace of Emitter @ SNR = %.2f',s2n(1));
plot(PL_Low,'r'); hold on; plot(BG_Low,'k'); 
ylim([0 1.5*max(PL_Low)]); legend('PL','BG'); title(tit);
subplot(2,1,2); tit = sprintf('Trace of Emitter @ SNR = %.2f', s2n(2));
plot(PL_High,'r'); hold on; plot(BG_High,'k'); 
ylim([0 1.5*max(PL_High)]);legend('PL','BG'); title(tit);

%% 2.7 blinks per second
offFrames2 = round(frames*rand(1,numOff(2))); % frames off @ 2.7bps
blinkLow2 = GPLow;  blinkHigh2 = GPHigh;

blinkLow2(:,:,offFrames2) = repmat(offLow,1,1,numOff(2));
    % 2bps @ SNR L
blinkHigh2(:,:,offFrames2) = repmat(offHigh,1,1,numOff(2));
    % 2bps @ SNR H

fLow2 = fullfile(F_save,'LowSNR_2o7bps.ome.tiff');
bfsave(uint16(blinkLow2),fLow2);
fHigh2 = fullfile(F_save,'HighSNR_2o7bps.ome.tiff');
bfsave(uint16(blinkHigh2),fHigh2);

%% Check 2.7bps trace
PL_Low2 = mean(reshape(blinkLow2(repmat(ROI_PL,1,1,frames)),A(1),frames));
PL_High2 = mean(reshape(blinkHigh2(repmat(ROI_PL,1,1,frames)),A(1),frames));
    % reshape = A-by-frames array where each column represents a frame
    % mean = 1-by-frames array where each column is the grey level of frame    
s2n = [mean(PL_Low2./BGSD_Low) mean(PL_High2./BGSD_High)]; % SNR of the traces
figure(n); n=n+1;
subplot(2,1,1); tit = sprintf('Trace of Emitter @ SNR = %.2f',s2n(1));
plot(PL_Low2,'r'); hold on; plot(BG_Low,'k'); 
ylim([0 1.5*max(PL_Low2)]); legend('PL','BG'); title(tit);
subplot(2,1,2); tit = sprintf('Trace of Emitter @ SNR = %.2f', s2n(2));
plot(PL_High2,'r'); hold on; plot(BG_High,'k'); 
ylim([0 1.5*max(PL_High2)]);legend('PL','BG'); title(tit);