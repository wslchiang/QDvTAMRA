function [F_save,BG,PL,cluster,skips] = getData(F,bits,fps,start,source,subBG)

% Generate general save file path based on date
formatOut = 'mmmm-dd-yyyy';
F_save = fullfile(F,datestr(date,formatOut));


% File path naming scheme for background and quantum dot PL
F_bg = fullfile(F,'BG-');
F_pl = fullfile(F,'PL-'); 
% Initializing cell arrays to store data of different sizes
BG = {}; PL = {}; sigma_BG = {}; sigma_PL = {};
cluster = []; skips = [];
fig1 = figure(1); fig2 = figure(2);

% Loop through data and do initial analysis w/ background subtraction
if(subBG)
    F_save = fullfile(F_save,'BG Subtracted'); %Update save directory
    % Create specific save directory and output files
    if(~exist(F_save,'dir'))
        mkdir(F_save);
    end
    filePL = fullfile(F_save,'PL.csv'); fidPL = fopen(filePL,'wt');
    fileBG = fullfile(F_save,'BG.csv'); fidBG = fopen(fileBG,'wt'); 
    fileBGSD = fullfile(F_save,'SD_BG.csv'); fidBGSD = fopen(fileBGSD,'wt');
    filePLSD = fullfile(F_save,'SD_PL.csv'); fidPLSD = fopen(filePLSD,'wt');
    fileClust = fullfile(F_save,'ClusterDetermination.csv'); 
    fidClust = fopen(fileClust,'wt');
    fprintf(fidClust, 'Cluster? (1=T), Low/Off Value, FWHM=2.355*std\n');
        
    for i = 1:length(dir(fullfile(F,'BG*.csv')))
        clf(fig1); clf(fig2); clf(figure(3));% clear figures at each iteration
        % Specific name for each individual file
        f_BG = [F_bg,num2str(i+start),'.csv'];
        f_PL = [F_pl,num2str(i+start),'.csv'];
        if isfile(f_BG) && isfile(f_PL) % Make sure files exist
            % Extract and save traces
            %tempBG = transpose(getBkg(f_BG));
            [tempBG,tempBG_SD] = getBkgSD(f_BG);
            tempBG = tempBG'; tempBG_SD = tempBG_SD';
            BG{i} = num2cell(tempBG,[1 2]);
            sigma_BG{i} = num2cell(tempBG_SD,[1 2]);
            % tempPL = transpose(getPL(f_PL));
            [tempPL,tempPL_SD] = getPLSD(f_PL);
            tempPL = tempPL'; tempPL_SD = tempPL_SD';
            PL{i} = num2cell(tempPL,[1 2]);
            sigma_PL{i} = num2cell(tempPL_SD,[1 2]);         
            
            %subtract BG
            noBG = tempPL-tempBG; noBG(noBG<0)=0;
            % Construct PDF from background subtracted PL
            if(bits == 8)
                [N,edges] = histcounts(noBG,'BinMethod','integers','Normalization','pdf');
            else
                [N,edges] = histcounts(noBG,'BinWidth',50,'Normalization','probability');
            end
            
            % Plot
            frames = length(noBG); x=linspace(0,frames/fps,frames);
            figure(1);
            subplot(1,2,1); plot(x,noBG); 
            ylim([0,2*max(noBG)]); legend('PL-BG');
            xlabel('Time(s)');ylabel('Grayscale Intensity');
            title('Background Subtracted PL');
            subplot(1,2,2); plot(N,edges(1:end-1)); hold on;
            
            %Fit PDF to determine single or cluster
            figure(2);
            h = histfit(noBG,length(edges),'kernel');
            figure(1); subplot(1,2,2); hold on;
            plot(h(2).YData/length(noBG),h(2).XData,'r');  
            legend('PDF','Fit');
            ylim([0,2*max(noBG)]); x_bound = maxk(N,2);
            xlim([min(N), min(x_bound)]);
            title('PDF of Grayscale Intensities');
            img = source+"-"+(i+start); 
            sgtitle(img); 
            name = fullfile(F_save,img+".png");
            exportgraphics(fig1,name);
            [pks,loc,wid] = findpeaks(h(2).YData,h(2).XData,'MinPeakProminence',0.125);            
            numPks = length(pks); 
            %figure(3); %findpeaks(h(2).YData,h(2).XData,'Annotate','extents');
            % [P,A] = autopeaks(h(2).XData,h(2).YData);
            %hits = findpeaks(h(2).YData,h(2).XData,'MinPeakProminence',0.125);
            %[fits,~,~,~,~,xi,yi] = peakfit([h(2).XData;h(2).YData],0,0,length(hits),3);
                % shapes (last input): 1=Gaussian, 2=Lorentzian, 3=Logistic
                % fits = pk#, position, height, width, area
                % can fit individual pks with xi and yi and findpeaks on
                % each to get widths if I want
            %name = fullfile(F_save,"PDF Fit of "+img+".png");
            %exportgraphics(figure(3),name);
            %numPks = size(fits,1);
            %loc = fits(:,2); pks = fits(:,3); wid = fits(:,4);
            if(numPks>2)
                % Multiple "peaks" -> Find and use largest 2
                [~,ind] = maxk(pks,2);
                if(sum(wid>(2^bits)/(10^(bits/8))==0))
                    % "narrow" peak defined at 10% of bit depth
                    cluster(i,1) = false;
                else
                    cluster(i,1) = true;
                end
                cluster(i,2) = loc(min(ind));
                cluster(i,3) = wid(min(ind));
            elseif(numPks==2)
                % 2 peaks
                [low,ind] = min(loc);
                if(sum(wid>(2^bits)/(10^(bits/8))==0))
                    % Peaks are distinct -> not cluster
                    cluster(i,1) = false;
                else
                    % Peaks are broad -> likely a cluster
                    cluster(i,1) = true;
                end
                cluster(i,2) = low;
                cluster(i,3) = wid(ind);
            elseif(numPks==1)                
                % LR = skewness(h(2).YData);
                    % negative if spreads left of mean
                    % right if spread right of mean
                    % NOT useful for my dataset
                LR = mean(noBG)-median(noBG);
                    % Positive skew (left peak) = mean > median
                if(wid<(2^bits)/(10^(bits/8)))
                    % narrow peak -> single QD
                    cluster(i,1) = false;
                else
                    % broad peak -> cluster
                    cluster(i,1) = true;
                end
                if(LR>0)
                    cluster(i,2) = loc;
                    cluster(i,3) = wid;
                else
                    cluster(i,2) = nan;
                    cluster(i,3) = nan;
                end
            else
                %return
                % unlikely to have a situation where no peak exists
                cluster(i,:) = repmat(nan,1,3);
                skips = [skips i];
            end
            % Export extracted values for future access
            fprintf(fidClust,'%d,%.2f,%.2f\n',cluster(i,:));
            fprintf(fidPL,repmat('%.2f,',1,length(tempPL)),tempPL); 
            fprintf(fidBG,repmat('%.2f,',1,length(tempBG)),tempBG);
            fprintf(fidBGSD,repmat('%.2f,',1,length(tempBG_SD)),tempBG_SD);
            fprintf(fidPLSD,repmat('%.2f,',1,length(tempPL_SD)),tempPL_SD);
            fprintf(fidPL,'\n'); fprintf(fidBG,'\n'); 
            fprintf(fidBGSD,'\n'); fprintf(fidPLSD,'\n');
        else
            skips = [skips i];
            BG{i} = {[]}; PL{i} = {[]}; 
            sigma_BG{i} = {[]}; sigma_PL{i} = {[]};
            cluster(i,:) = repmat(nan,1,3);
            fprintf(fidPL,'%.2f\n',[]); fprintf(fidBG,'%.2f\n',[]);
            fprintf(fidBGSD,'%.2f\n',[]); fprintf(fidPLSD,'%.2f\n',[]);
            fprintf(fidClust,'%d,%.2f,%.2f\n',cluster(i,:));
        end
    end
else % Same set of commands, except no BG subtract
    % Specific save location
    F_save = fullfile(F_save,'Raw Data');
    if(~exist(F_save,'dir'))
        mkdir(F_save); 
    end
    % Output files
    filePL = fullfile(F_save,'PL.csv'); fidPL = fopen(filePL,'wt');
    fileBG = fullfile(F_save,'BG.csv'); fidBG = fopen(fileBG,'wt'); 
    fileBGSD = fullfile(F_save,'SD_BG.csv'); fidBGSD = fopen(fileBGSD,'wt');
    fileClust = fullfile(F_save,'ClusterDetermination.csv'); 
    fidClust = fopen(fileClust,'wt');
    fprintf(fidClust, 'Cluster? (1=T), Low/Off Value, FWHM=2.355*std\n');
    
    for i = 1:length(dir(fullfile(F,'BG*.csv')))
        clf(fig1); clf(fig2); % clear figure
        % Specific name for each individual file
        f_BG = [F_bg,num2str(i+start),'.csv'];
        f_PL = [F_pl,num2str(i+start),'.csv'];
        if isfile(f_BG) && isfile(f_PL) %Make sure files exist
            % Get and store BG & PL traces
            [tempBG,tempBG_SD] = getBkgSD(f_BG);
            [tempPL,tempPL_SD] = getPLSD(f_PL);
            tempBG = tempBG'; tempBG_SD = tempBG_SD';
            tempPL = tempPL'; tempPL_SD = tempPL_SD';
            BG{i} = num2cell(tempBG,[1 2]);
            PL{i} = num2cell(tempPL,[1 2]);
            sigma_BG{i} = num2cell(tempBG_SD,[1 2]);
            sigma_PL{i} = num2cell(tempPL_SD,[1 2]);
            % Construct PDF for PL
            [N,edges] = histcounts(tempPL,'BinMethod','integers','Normalization','pdf');
            % Plot
            frames = length(tempPL); x = linspace(0,frames/fps,frames);
            figure(1);
            subplot(1,2,1); plot(x,tempPL); hold on; plot(x,tempBG);
            xlabel('Time (s)'); ylabel('Grayscale Intensity');
            title('Raw PL & BG'); legend('PL','BG'); ylim([0,2^bits-1]);
            subplot(1,2,2); plot(N,edges(1:end-1)); hold on;
                       
            % determine if PDF is
            % bimodal (single)
            % or broad (cluster)
            figure(2);
            h = histfit(tempPL,length(edges),'kernel');
            figure(1); subplot(1,2,2); hold on;
            plot(h(2).YData/length(tempPL),h(2).XData,'r');  
            legend('PDF','Fit'); 
            ylim([0,2*max(tempPL)]); x_bound = maxk(N,2);
            xlim([min(N), min(x_bound)]);
            % ylim([0,2^bits-1]);
            title('PDF of Grayscale Intensities');
            img = source+"-"+(i+start); sgtitle(img); name = fullfile(F_save,img+".png");
            exportgraphics(fig1,name);
            [pks,loc,wid] = findpeaks(h(2).YData,h(2).XData,'MinPeakProminence',0.125);
            numPks = length(pks); figure(3);
            findpeaks(h(2).YData,h(2).XData,'Annotate','extents');
            %[P,A] = autofindpeaks(h(2).XData,h(2).YData);
            %[fits,~,~,~,~,xi,yi] = peakfit([h(2).XData;h(2).YData],0,0,size(P,1),2);
                % fits = pk#, position, height, width, area
                % can fit individual pks with xi and yi and findpeaks on
                % each to get widths if I want
            %numPks = size(fits,1);
            %loc = fits(:,2); pks = fits(:,3); wid = fits(:,4);
            if(numPks>2)
                % Multiple "peaks" -> Find and use largest 2
                [~,ind] = maxk(pks,2);
                if(sum(wid>(2^bits)/(10^(bits/8)))==0)
                    % "narrow" peak defined at 10% of bit depth
                    cluster(i,1) = false;
                else
                    cluster(i,1) = true;
                end
                cluster(i,2) = loc(min(ind));
                cluster(i,3) = wid(min(ind));
            elseif(numPks==2)
                % 2 peaks
                [low,ind] = min(loc);
                if(sum(wid>(2^bits)/(10^(bits/8)))==0)
                    % Peaks are distinct -> not cluster
                    cluster(i,1) = false;
                else
                    % Peaks are broad -> likely a cluster
                    cluster(i,1) = true;
                end
                cluster(i,2) = low;
                cluster(i,3) = wid(ind);
            elseif(numPks==1)                
                % LR = skewness(h(2).YData);
                    % negative if spreads left of mean
                    % right if spread right of mean
                    % NOT useful for my dataset
                LR = mean(tempPL)-median(tempPL);
                    % Positive skew (left peak) = mean > median
                if(wid<(2^bits)/(10^(bits/8)))
                    % narrow peak -> single QD
                    cluster(i,1) = false;
                else
                    % broad peak -> cluster
                    cluster(i,1) = true;
                end
                if(LR>0)
                    cluster(i,2) = loc;
                    cluster(i,3) = wid;
                else
                    cluster(i,2) = nan;
                    cluster(i,3) = nan;
                end
            else
                return
                % unlikely to have a situation where no peak exists
            end
            % Output extracted data to a file to access later
            fprintf(fidClust,'%d,%.2f,%.2f\n',cluster(i,:));
            fprintf(fidPL,repmat('%.2f,',1,length(tempPL)),tempPL); 
            fprintf(fidBG,repmat('%.2f,',1,length(tempBG)),tempBG);
            fprintf(fidBGSD,repmat('%.2f,',1,length(tempBG_SD)),tempBG_SD);
            fprintf(fidPLSD,repmat('%.2f,',1,length(tempPL_SD)),tempPL_SD);
            fprintf(fidPL,'\n'); fprintf(fidBG,'\n'); fprintf(fidBGSD,'\n');
        else
            skips = [skips i];
            BG{i} = {[]}; PL{i} = {[]}; 
            sigma_BG{i} = {[]}; sigma_PL{i} = {[]};
            cluster(i,:) = repmat(nan,1,3);
            fprintf(fidPL,'%.2f\n',[]); fprintf(fidBG,'%.2f\n',[]);
            fprintf(fidBGSD,'%.2f\n',[]); fprintf(fidPLSD,'%.2f\n',[]);
            fprintf(fidClust,'%d,%.2f,%.2f\n',cluster(i,:));
        end
    end
end
fclose('all');

