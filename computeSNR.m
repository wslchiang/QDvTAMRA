%% compareSNR
% v2.0 | compare mean PL to std BG
% August 23, 2020
%% Select File Locations
F_QD = uigetdir(pwd,'Select Folder for QD Files');
F_TAMRA = uigetdir(pwd,'Select Folder for TAMRA Files');
    % File location/folder where the extracted and compiled BG and PL csv
    % files are
%% Extract values
F_QD_PL = fullfile(F_QD,'PL.csv');
F_QD_SDBG = fullfile(F_QD,'SD_BG.csv');
F_QD_SDPL = fullfile(F_QD,'SD_PL.csv'); 
F_TAMRA_PL = fullfile(F_TAMRA,'PL.csv');
F_TAMRA_SDBG = fullfile(F_TAMRA,'SD_BG.csv');
F_TAMRA_SDPL = fullfile(F_TAMRA,'SD_PL.csv');

if(isfile(F_QD_PL) && isfile(F_QD_SDBG) && isfile(F_QD_SDPL))
    PL_QD = readmatrix(F_QD_PL); PL_QD = PL_QD(:,1:end-1);
    SDBG_QD = readmatrix(F_QD_SDBG); SDBG_QD = SDBG_QD(:,1:end-1);
    SDPL_QD = readmatrix(F_QD_SDPL); SDPL_QD = SDPL_QD(:,1:end-1);
else
    error('Necessary QD files do not exist in this directory \n');
    return
end

if(isfile(F_TAMRA_PL) && isfile(F_TAMRA_SDBG) && isfile(F_TAMRA_SDPL))
    PL_TAMRA = readmatrix(F_TAMRA_PL); PL_TAMRA = PL_TAMRA(:,1:end-1);
    SDBG_TAMRA = readmatrix(F_TAMRA_SDBG); 
    SDBG_TAMRA = SDBG_TAMRA(:,1:end-1);
    SDPL_TAMRA = readmatrix(F_TAMRA_SDPL);
    SDPL_TAMRA = SDPL_TAMRA(:,1:end-1);
else
    error('Necessary TAMRA files do not exist in this directory \n');
    return
end
%%
fileQD = fullfile(F_QD,'SNR2.csv');
fidQD = fopen(fileQD,'wt');
fileTAMRA = fullfile(F_TAMRA,'SNR2.csv');
fidTAMRA = fopen(fileTAMRA,'wt');
%% Compute individual SNR for QD
R_QD = repmat(0,1,size(PL_QD,1));
    % Initialize storage for SNR values with zeroes to decrease overhead
for i=1:size(PL_QD,1)
    %tempPL = PL_QD(i,:); tempSD = SD_QD(i,:);
    %R_QD(i) = mean(tempPL)/mean(tempSD);
    R_QD(i) = mean(PL_QD(i,:)./SDBG_QD(i,:));
    fprintf(fidQD,'%.2f\n',R_QD(i));    
end
%% Statistical Analysis and plotting for QD
[N,edges] = histcounts(R_QD,'BinMethod','integers');
centers = mean([edges(1:end-1);edges(2:end)]);
mu_QD = mean(R_QD); sigma_QD = std(R_QD);
fprintf(fidQD,'average, %.2f\n', mu_QD);
fprintf(fidQD,'std, %.2f\n', sigma_QD);
fig = figure;
bar(centers,N); ylabel('Occurences'); xlabel('SNR'); 
tit = sprintf('SNR of QDs \\mu = %.2f \\sigma = %.2f', [mu_QD sigma_QD]);
title(tit);
name = fullfile(F_QD,'SNR2 of QDs.png');
exportgraphics(fig,name);
%% Compute individual SNR for TAMRA
R_TAMRA = repmat(0,1,size(PL_TAMRA,1));
    % Initialize storage for SNR values with zeroes to decrease overhead
for i=1:size(PL_TAMRA,1)
    % R_TAMRA(i) = mean(PL_TAMRA(i,:))/mean(SD_TAMRA(i,:));
    R_TAMRA(i) = mean(PL_TAMRA(i,:)./SDBG_TAMRA(i,:));
    fprintf(fidTAMRA,'%.2f\n',R_TAMRA(i));
end
%% Statistical Analysis and plotting for TAMRA
[N,edges] = histcounts(R_TAMRA,'BinMethod','integers');
centers = mean([edges(1:end-1);edges(2:end)]);
mu_TAMRA = mean(R_TAMRA); sigma_TAMRA = std(R_TAMRA);
fprintf(fidTAMRA,'average, %.2f\n', mu_TAMRA);
fprintf(fidTAMRA,'std, %.2f\n', sigma_TAMRA);
fig2 = figure;
bar(centers,N); ylabel('Occurences'); xlabel('SNR'); 
tit = sprintf('SNR of TAMRAs \\mu = %.2f \\sigma = %.2f', [mu_TAMRA sigma_TAMRA]);
title(tit);
name2 = fullfile(F_TAMRA,'SNR2 of TAMRAs.png');
exportgraphics(fig2,name2);
%% End
fprintf('SNR of QD = %.2f +/- %.2f\n',[mu_QD sigma_QD]);
fprintf('SNR of TAMRA = %.2f +/- %.2f\n',[mu_TAMRA sigma_TAMRA]);
fclose('all')
close all
%clear all