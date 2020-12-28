%{
    v1.2 | october 26, 2020
    Created by Wesley Chiang
    Modified for analysis of multiemitter outputs
%}
%% Define Parameters
clear all
temp = load('truePx.mat');
true_px = temp.true_px;
    % true positions of the emitter in px
conversion = input('What is the pixel to nm conversion? ');
    % user input for calculated pixel size based on system
true_nm = true_px*conversion;   
    % True positions in nm
%% Create Save File
path = uigetdir(pwd);
F_save = fullfile(path,'Results.csv'); fid = fopen(F_save,'wt');
head = 'SNR,bps,QuickPALM,Cluster,MuX (px),MuY (px),StdX (px),StdY (px)\n';
fprintf(fid,head);
%% Loop
more = true;
while(more)
    fignum = 1;
    %% File Selection & Extract Data
    [file,folder] = uigetfile(fullfile(path,'*.csv'),'Select csv file to analyze');
    s2n = folder(end-12:end-9); bps = folder(end-4);
    if(strcmpi(s2n(1),'\'))
        s2n=s2n(2:end);
    end
    f = fullfile(folder,file);
    [Xpx,Ypx,Xnm,Ynm,frame] = getData(f);
    spots = [Xpx Ypx]+1;

    %% Define Analysis
    opt = input('Which Cluster Method? (I)nconsistent or (M)ax Cluster? ','s');
    %% Cluster Analysis
    c = getClust(spots,opt,fignum,f(1:end-4));
    idx = str2num(input('Enter space separated values clusters to use ','s'));
        % 1 = bottom left, 2 = middle right, 3 = top left
    for i=1:length(idx)
        vals = spots(c==idx(i),:);
        if(isempty(vals))
            fprintf('No points at this cluster position \n')
            output = [s2n,',',bps,',',file,',',opt,',%.3f,%.3f,%.3f,%.3f\n'];
            fprintf(fid,output,[nan,nan,nan,nan])
        else
            [muX,muY,sigmaX,sigmaY] = CoM(vals,opt,true_px,i,fignum,f(1:end-4));
            output = [s2n,',',bps,',',file,',',opt,',%.3f,%.3f,%.3f,%.3f\n'];
            fprintf(fid,output,[muX,muY,sigmaX,sigmaY])
        end
    end
    more = strcmpi(input('Analyze more files (Y/N)? ','s'),'y');
end
%% End
fclose('all')

%% Get Clusters
    % Adapted from MultiEmitter2.m
function c = getClust(spots,opt,fignum,fp)
    D = pdist(spots);
    L = linkage(D,'average');
    coeff = cophenet(L,D)
    I = inconsistent(L);
    if(strcmpi(opt,'I'))
        cutoff = mean(I(:,4))+1.25*std(I(:,4))
        c = clusterdata(spots,'cutoff',cutoff);
        fig = figure(fignum); 
        clf
        fignum=fignum+1;
        gscatter(spots(:,1),spots(:,2),c,parula(max(c)));
        title('Cluster Analysis Using Inconsistency Cutoff')
    elseif(strcmpi(opt,'M'))
        c = cluster(L,'maxclust',7);
        fig = figure(fignum); 
        clf
        fignum=fignum+1; 
        gscatter(spots(:,1),spots(:,2),c,parula(max(c)));
        title('Cluster Analysis Using MaxClust=7')
    else
        fprintf('Entered Invalid Option for Cluster Analysis\n')
    end
    fprintf(['Cophenetic Correlation of cluster = ',num2str(coeff),'\n'])
    figname = [fp,'_',opt,'.png'];
    exportgraphics(fig,figname);
end
%% Fit Gaussians
function [muX,muY,sigmaX,sigmaY] = CoM(vals,opt,true_px,i,fignum,fp)
    if(size(vals,1)==1)
        muX = vals(1,1); muY = vals(1,2);
        sigmaX = nan; sigmaY = nan;
    else
        [~,edgeX]=histcounts(vals(:,1)); [~,edgeY]=histcounts(vals(:,2));
        fig = figure(fignum);
        clf
        fignum=fignum+1; 
        subplot(1,2,1); hx=histfit(vals(:,1),length(edgeX));
        xlabel('X (px)'); ylabel('Counts');
        [Xpks,Xloc,Xwid] = findpeaks(hx(2).YData,hx(2).XData);
        muX = Xloc; sigmaX = Xwid/2.35;
        leg = sprintf('\\mu = %.2f +/- \\sigma = %.2f',[muX,sigmaX]);
        legend(hx(2),leg)
        subplot(1,2,2); hy=histfit(vals(:,2),length(edgeY));
        xlabel('Y (px)');
        [Ypks,Yloc,Ywid] = findpeaks(hy(2).YData,hy(2).XData);
        muY = Yloc; sigmaY = Ywid/2.35;
        leg = sprintf('\\mu = %.2f +/- \\sigma = %.2f',[muY,sigmaY]);
        legend(hy(2),leg)
        true_loc = true_px(i,:);
        tit = sprintf('Gaussian Fit for Centroid at Position (%.2f,%.2f)',true_loc);
        sgtitle(tit);
        figname = [fp,'_'+opt,'_Cluster',num2str(i),'.png'];
        exportgraphics(fig,figname);
    end
end

    