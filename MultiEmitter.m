%% Simulate Emitters
% v3.1 | Multiple emitters in ROI
% based on SNR = mean(PL)/std(BG) (v2.0)
% Introduce more noise/variance between frames (v2.1)
% September 8, 2020
% by Wesley Chiang
% Target SNR QD=5.33 | R6G = 2.02
%% Parameter setup
% Single emitter parameters
ROI = 8;        % 8x8 pixel ROI equivalent to that used to get traces
ROI_PL = repmat(false,ROI,ROI); 
ROI_PL(2:end-1,2:end-1)=true;  % 6x6 region used to get PL traces
ROI_BG = ~ROI_PL;           % region wrapping around PL ROI to get BG trace
% Multiemitter parameters
ROI_2 = 11;     % 11x11 pixel ROI to create overlapped emitters
startX = [2, 4, 1];         % TL x-position of overlapped emitters (col)
startY = [1, 3, 4];         % TL y-position of overlapped emitters (row)
    % MATLAB indices go row,col
A = [sum(sum(ROI_PL)) sum(sum(ROI_BG))];    % A represents are of ROI
[xPL,yPL] = find(ROI_PL);   % coordinates for PL ROI
[xBG,yBG] = find(ROI_BG);   % coordinates for BG ROI
bits = 16;      % bit depth used for image
QE = 0.8;        % Quantum effiency of pco edge 4.2 @ 600nm = 80%
depth = 45000;  % Quantum Well Depth of pco edge 4.2 = 45000 e-
noise = 2.3;    % rms readout noise of pco edge 4.2 = 2.3 e-
fps = 40;       % frame rate of traces

%{
%% Blink Parameters
bps = [2,2.7];                  % blinks per second
numOff = round(bps/40*frames);  % number of frames to be off
%}

%% Save Location
F_save = uigetdir(pwd,'Set Save Location');

%% Create emitter of desired SNR
HoL = input('High or Low SNR? ','s');
m = [1:ROI]-(ROI+1)/2;  % grid layout of length ROI centered at 0
[x,y] = meshgrid(m,m);  % m-by-m grid of values
ind = x.^2+y.^2;        % square of indices of size ROI
emit = ind<((ROI/4)^2);
ap = ind<((ROI/8)^2);
n = 1; 
%{
figure(n); n=n+1; 
subplot(2,1,1); imagesc(emit); colormap(gray);axis image; title('Emitter');
subplot(2,1,2); imagesc(ap); colormap(gray); axis image; title('Aperture');
%}
clim = [0 15000];       % Scale upper and lower bounds for plot

if(strcmpi(HoL,'Low'))
    emitter = 500*double(emit); % emitter low SNR
    emitter(~emit)=1; %emit8 = 500*emit8;
    %{
    figure(n); n=n+1; 
    imagesc(emitter); title('Emitter Low SNR');
    colormap(gray); colorbar; axis image
    %}
elseif(strcmpi(HoL,'High'))
    emitter = 350*double(emit);% emitter high SNR
    emitter(~emit)=1; %emit3 = 500*emit3;
    %{
    figure(n); n =n+1;
    imagesc(emitter); title('Emitter High SNR');
    colormap(gray); colorbar; axis image
    %}
else
    return
end
%% Add blur to emitters
radcoeff = 2*pi/((ROI/16)^2); q_phase = exp(1i*radcoeff*ind);
    % Generate quadrat
h = fftshift(fftn(ifftshift(ap.*q_phase)));
PSF = abs(h).^2;
H = fftshift(fftn(ifftshift(PSF)));

F = fftshift(fftn(ifftshift(emitter)));
B = F.*H;
b = fftshift(ifftn(ifftshift(B)));
r = mean(b(ROI_PL))/std(b(ROI_BG));

figure(n); n=n+1;
imagesc(1/4*b);
r = mean(b(ROI_PL))/std(b(ROI_BG));
tit = sprintf('Blurred emitter SNR = %.2f',r);
title(tit);
colormap(gray); colorbar; axis image

%% Create n-frame movies of emitters
frames = input('Number of frames? ');  % number of frames used in traces
mov = repmat(1/3*b,1,1,frames);
%% Setup Multi-emitter movie
numEmit = input('Number of Emitters in this simulation? ');
bps = str2num(input('Enter space separated values for blink rate ','s'));
numOff = round(bps/40*frames);  % number of frames to be off
if numEmit > 1
    movROI = repmat(zeros(ROI_2),1,1,frames);
else
    movROI = repmat(zeros(ROI),1,1,frames);
end
% Introduce desired SNR

for i = 1:length(bps)
    full{i} = movROI;
end

   
if(strcmpi(HoL,'Low'))
    scaleBG = 20/r*(-0.5+(1.25+0.5)*rand(size(mov)));
    mov_s2n = mov.*ROI_BG.*scaleBG+mov; mov_s2n(mov_s2n<0) = 0;
        % Increase noise in BG to achieve desired SNR
    PL = mean(reshape(mov_s2n(repmat(ROI_PL,1,1,frames)),A(1),frames));
    BG = mean(reshape(mov_s2n(repmat(ROI_BG,1,1,frames)),A(2),frames));
    BGSD = std(reshape(mov_s2n(repmat(ROI_BG,1,1,frames)),A(2),frames));
    s2n = mean(PL./BGSD);
    scalePL = sqrt(s2n)*(0.35+(1-0.35)*rand(size(mov_s2n)));
    mov_s2n = 2/3*mov_s2n.*ROI_PL.*scalePL+mov_s2n.*ROI_BG;
        % Introduce random noise across emitter pixels to create frame to frame
        % variance that ~SNR
    PL = mean(reshape(mov_s2n(repmat(ROI_PL,1,1,frames)),A(1),frames));
    s2n = mean(PL./BGSD);
elseif(strcmpi(HoL,'High'))
    % High SNR
    scaleBG = 4/r*(-1.5+(1.75+1.5)*rand(size(mov)));
    mov_s2n = mov.*ROI_BG.*scaleBG+mov;
        % Noise in BG
    PL = mean(reshape(mov_s2n(repmat(ROI_PL,1,1,frames)),A(1),frames));
    BG = mean(reshape(mov_s2n(repmat(ROI_BG,1,1,frames)),A(2),frames));
    BGSD = std(reshape(mov_s2n(repmat(ROI_BG,1,1,frames)),A(2),frames));
    s2n = mean(PL./BGSD);
    scalePL = 1/2*sqrt(s2n)*(0.4+(1-0.4)*rand(size(mov_s2n)));
    mov_s2n = mov_s2n.*ROI_PL.*scalePL+mov_s2n.*ROI_BG;
        % Frame to frame variance
    PL = mean(reshape(mov_s2n(repmat(ROI_PL,1,1,frames)),A(1),frames));
    s2n = mean(PL./BGSD);
end

% Introduce Poisson-Gaussian Noise due to camera parameters
Poiss = poissrnd(round(mov_s2n*QE*depth/(2^bits)));
    % Poisson (shot noise)
    % Random poisson noising ontop of photoelectrons from simulated emitter
    % Quantum efficiency and well depth of pco edge 4.2
PoissGauss = Poiss+sqrt(noise)*randn(size(Poiss)); 
    % Gaussian + Poisson noise

blinking = PoissGauss;
BG = mean(reshape(blinking(repmat(ROI_BG,1,1,frames)),A(2),frames));
BGSD = std(reshape(blinking(repmat(ROI_BG,1,1,frames)),A(2),frames));
    % reshape = A-by-frames array where each column represents a frame
    % mean = 1-by-frames array where each column is the grey level of frame
off = mean(BG)+sqrt(mean(BGSD))*randn(size(ROI_PL)); 
off(off<0)=0;

%% need to add loop to add 3 of the same emitters
for i = 1:numEmit
    for b = 1:length(numOff)
        blinking = PoissGauss;
        offFrames = round(frames*rand(1,numOff(b))); % frames off 
        blinking(:,:,offFrames) = repmat(off,1,1,numOff(b));
            % set randomly selected frames to be off for bps(i)
            % Set coordinates of PL ROI to the average BG level @ SNR  
        if numEmit == 1
            full{b} = blinking;
        else
            cols = startX(i):startX(i)+ROI-1; % x vals
            rows = startY(i):startY(i)+ROI-1; % y vals
            temp = full{b}(rows,cols,:);
            full{b}(rows,cols,:)=temp+blinking;
        end
    end
end
for b=1:length(bps)
    fileName = HoL+"SNR_"+bps(b)+"bps.ome.tiff";
    f = fullfile(F_save,char(fileName));
    bfsave(uint16(full{b}),f);
end


