%{
    v1.0 | September 2, 2020
    Created by Wesley Chiang
%}
%% Define Parameters
clear all
True_px = [4.5 , 4.5];
    % true positions of the emitter in px
conversion = input('What is the pixel to nm conversion? ');
    % user input for calculated pixel size based on system
True_nm = True_px*conversion;   
    % True positions in nm
F = uigetdir(pwd,'Select Folder w/ Subfolders');
rounds = input('How Many Rounds (Subfolders)? ');
F_save = fullfile(F,'Results.csv'); fid = fopen(F_save,'wt');
fprintf(fid,'SNR,bps,Mu_X (px),Mu_Y (px),Std_X (px),Std_Y (px),%%Err X,%%Err Y\n');

%% High SNR @ 2 blinks per second (bps)
for i = 1:rounds
    temp = "Round"+i;
    F_temp = fullfile(F,temp,'HighSNR_2bps.csv');
        % File path for High SNR 2bps
    [Xpx, Ypx, Xnm, Ynm, FrameNumber]=getData(F_temp);
        % Import Data
    [muX,muY,sigmaX,sigmaY] = CoM(Xpx,Ypx);
        % Fit ~ Gaussian parameters
    mu_px = [muX,muY]+1; sigma_px = [sigmaX,sigmaY];
        %Account for ImageJ indexing @ 0 vs MATLAB @ 1
    err_px = (mu_px-True_px)./True_px*100;   
        % Percent error in simulation
    fprintf(fid,'High,2,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',[mu_px,sigma_px,err_px]);
end
%% Low SNR @ 2bps
for i=1:rounds
    temp = "Round"+i;
    F_temp = fullfile(F,temp,'LowSNR_2bps.csv');
        % File path for Low SNR 2bps
    [Xpx, Ypx, Xnm, Ynm, FrameNumber]=getData(F_temp);
        % Import Data
    [muX,muY,sigmaX,sigmaY] = CoM(Xpx,Ypx);
        % Fit ~ Gaussian parameters
    mu_px = [muX,muY]+1; sigma_px = [sigmaX,sigmaY];
        % Account for ImageJ indexing @ 0 vs MATLAB @ 1
    err_px = (mu_px-True_px)./True_px*100;   
        % Percent error in simulation
    fprintf(fid,'Low,2,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',[mu_px,sigma_px,err_px]);
end
%% High SNR @ 2.7bps
for i=1:rounds
    temp = "Round"+i;
    F_temp = fullfile(F,temp,'HighSNR_2o7bps.csv');
        % File path for High SNR 2.7bps
    [Xpx, Ypx, Xnm, Ynm, FrameNumber]=getData(F_temp);
        % Import Data
    [muX,muY,sigmaX,sigmaY] = CoM(Xpx,Ypx);
        % Fit ~ Gaussian parameters
    mu_px = [muX,muY]+1; sigma_px = [sigmaX,sigmaY];
        % Account for ImageJ indexing @ 0 vs MATLAB @ 1
    err_px = (mu_px-True_px)./True_px*100;   
        % Percent error in simulation
    fprintf(fid,'High,2.7,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',[mu_px,sigma_px,err_px]);
end
%% Low SNR @ 2.7bps
for i=1:rounds
    temp = "Round"+i;
    F_temp = fullfile(F,temp,'LowSNR_2o7bps.csv');
        % File path for Low SNR 2.7bps
    [Xpx, Ypx, Xnm, Ynm, FrameNumber]=getData(F_temp);
        % Import Data
   [muX,muY,sigmaX,sigmaY] = CoM(Xpx,Ypx);
        % Fit ~ Gaussian parameters
    mu_px = [muX,muY]+1; sigma_px = [sigmaX,sigmaY];
        % Account for ImageJ indexing @ 0 vs MATLAB @ 1
    err_px = (mu_px-True_px)./True_px*100;   
        % Percent error in simulation
    fprintf(fid,'Low,2.7,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',[mu_px,sigma_px,err_px]);
end
%% End
fclose('all')

%% Fit Gaussians
function [muX,muY,sigmaX,sigmaY] = CoM(Xpx,Ypx)

%{
---- for CoM ---
[N,edgeY,edgeX] = histcounts2(Xpx,Ypx);
    % bin detected pixel centers
    % input has X and Y swapped becauase of array indexing in histcounts2
centerX = mean([edgeX(1:end-1);edgeX(2:end)]); 
centerY = mean([edgeY(1:end-1);edgeY(2:end)]);
    % bin centers to match dimensions of hist matrix
Xcol = sum(N); Yrow = sum(N,2); Yrow = Yrow';
    % sum down columns for weighted X
    % sum across rows for weighted Y
total = sum(N,'all');
Xws = sum(Xcol.*centerX); Yws = sum(Yrow.*centerY);
    % Weighted sum to get center of mass
muX = Xws/total; muY = Yws/total;
    % center of mass in X and Y
sigmaX = sqrt(sum(Xcol.*abs(centerX-muX))/total);
sigmaY = sqrt(sum(Yrow.*abs(centerY-muY))/total);
    % standard deviation in localization
%}
Xpx = Xpx+1; Ypx = Ypx+1;
[~,edgeX]=histcounts(Xpx); [~,edgeY]=histcounts(Ypx);
figure; 
subplot(1,2,1); hx=histfit(Xpx,length(edgeX));
xlabel('X (px)'); ylabel('Counts');
[Xpks,Xloc,Xwid] = findpeaks(hx(2).YData,hx(2).XData);
subplot(1,2,2); hy=histfit(Ypx,length(edgeY));
xlabel('Y (px)');
[Ypks,Yloc,Ywid] = findpeaks(hy(2).YData,hy(2).XData);
muX = Xloc; sigmaX = Xwid/2.35;
muY = Yloc; sigmaY = Ywid/2.35;
sgtitle('Fitted Histograms for Centroid Location');
end
    