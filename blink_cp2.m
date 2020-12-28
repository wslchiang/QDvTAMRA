% Backgrouund subtraction, then changepoint
function [results,OnTime,OffTime]=blink_cp2(F_save,bits,fps,source,BG,PL,cluster,skips)
%Results array for blink, slope on, error, slope off, error
%Cell Arrays to store information about lifetimes 
    %col 1 = times, col 2 = prob
OnTime = {}; OffTime = {};
%Temporary on off arrays to store count different lifetimes
on = []; off = []; all_on = []; all_off = [];
s_on=[]; s_off=[]; % single QDs
clust_on=[];clust_off=[]; %clusters
%Temporary indexing variables
index = 0; events = 0; fluct = 0;
counts = [0 0]; %counts(1) = on, counts(2) = off
figure(3); figure(4); thresh = '';

F_save = fullfile(F_save, 'Changepoint+Blink');
if(~exist(F_save,'dir'))
    mkdir(F_save);
end
fileName = fullfile(F_save,'Results.csv');
fileName2 = fullfile(F_save,'Log.txt');
fid = fopen(fileName,'wt'); fid_log = fopen(fileName2,'wt');
fprintf(fid, 'Blinks per second, ON Fit Value, ON Fit Error, OFF Fit Value, OFF Fit Error\n');
for j = 1:length(BG)
    thresh ='';
    if any(skips(:) == j)
         fprintf(fid_log,"Skipped Data Set at interation " + j)
         results(j,:) = repmat(nan,1,5);
         fprintf(fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n',results(j,:));
         OnTime{j} = {[]}; OffTime{j} = {[]};
        continue
    end
    temp_PL = [PL{j}{:}];
    temp_BG = [BG{j}{:}];
    noBG = temp_PL-temp_BG; noBG(noBG<0)=0;
    %if(strcmpi(source,"QD"))
        if(cluster(j,1))
            if(sum(isnan(cluster(j,:)))==0)
                cutoff = cluster(j,2)+cluster(j,3)/2.355;
                thresh = 'Peak_{Low,PDF}+\sigma_{PDF}';
            else
                cutoff = min(temp_BG)+3*std(temp_BG);
                thresh = 'Min_{BG}+3\sigma_{BG}';
            end
        else
            cutoff = cluster(j,2)+3/2.355*cluster(j,3);
            thresh = 'Peak_{Low,PDF}+3\sigma_{PDF}';
        end
    %else
    %    cutoff = 3*std(temp_BG);
    %    thresh = '3\sigma_{BG}';
    %end
    [~,s] = ischange(noBG,'Threshold',2^(bits*1.25));
    OnOff = s >= cutoff;
    index = 1;
    on = []; off =[];
    while index < length(OnOff)-1
        counts = [0 0];
        while OnOff(index) == 1 && index ~= length(OnOff)
            counts(1) = counts(1)+1;
            index = index+1;
        end
        on = [on counts(1)];
        while OnOff(index) == 0 && index ~= length(OnOff)
            counts(2) = counts(2)+1;
            index = index+1;
        end
        off = [off counts(2)];
    end
    %Remove any potential leading zero values
    on = on(on~=0);
    off = off(off~=0);
    %Number of transitions/fluctuations in trace 
    fluct = (length(on)+length(off))/(length(OnOff)/fps); % Fluctuations per second
    %Bin different lifetimes where index = time on
    t_on = histcounts(on,'BinMethod','integers'); 
    t_off = histcounts(off,'BinMethod','integers');
    %Store indices of relevant on and off times
    on_val = find(t_on~=0)/fps; off_val = find(t_off~=0)/fps;
    if or(length(on_val)<2,length(off_val)<2)
        fprintf(fid_log,"Pair of data at " + j + " not usable\n");
        results(j,:) = repmat(nan,1,5);
        fprintf(fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n',results(j,:));
        OnTime{j} = {[]}; OffTime{j} = {[]};
        continue
    end
    events = length(on_val)+length(off_val); % number of unique lifetimes
    %remove nonoccuring on and off times
    t_on = t_on(t_on~=0); t_off = t_off(t_off~=0);
    %probability of on and off events
    p_on = t_on./events; p_off = t_off./events;
    %Probability Density of on and off events
    pdf_on = p_on./on_val; pdf_off = p_off./off_val;
    %Outputting essential values for further analysis
    if(cluster(j,1))
        clust_on = [clust_on [on_val;pdf_on]];
        clust_off = [clust_off [off_val;pdf_off]];
    else
        s_on = [s_on [on_val;pdf_on]];
        s_off = [s_off [off_val;pdf_off]];
    end
    % all_on = [all_on [on_val;pdf_on]];
    % all_off = [all_off [off_val;pdf_off]];
    
    OnTime{j} = [on_val' p_on'];
    OffTime{j} = [off_val' p_off'];
    %Fitting for power law
    [c_on,gof_on] = fit(log10(on_val'),log10(pdf_on'),'poly1');
    [c_off,gof_off] = fit(log10(off_val'),log10(pdf_off'),'poly1');
    m_on = sprintf('m_{on} = %.2f (%.2f)',[c_on.p1,gof_on.rmse]);
    m_off = sprintf('m_{off} = %.2f (%.2f)',[c_off.p1,gof_off.rmse]);
    %Plot the fit in log-log
    figure(3);
    fig = tiledlayout(1,2);
    nexttile;
    loglog(on_val,10.^polyval([c_on.p1,c_on.p2],log10(on_val)),'b');
    hold on
    loglog(on_val,pdf_on,'ok');
    legend(m_on); title('On Time'); 
    ylabel('log(PD) (1/s)'); xlabel('log(time) (s)');
    nexttile;
    loglog(off_val,10.^polyval([c_off.p1,c_off.p2],log10(off_val)),'b');
    hold on
    loglog(off_val,pdf_off,'ok');
    legend(m_off); title('Off Time'); xlabel('log(time) (s)');
    %Save the figure as png
    img = "Fit of "+ source + j + ".png"; name = fullfile(F_save,img);
    exportgraphics(fig,name);
    %Plot Original Data
    figure(4);
    img = "Blinking of " + source + j + ".png"; name = fullfile(F_save,img);
    subplot(2,1,1);
    hold on
    plot(linspace(0,length(noBG)/fps,length(noBG)),noBG,'b');
    stairs(linspace(0,length(noBG)/fps,length(noBG)),s,'r');
    yline(cutoff,'-.k'); ylim([0,2*max(noBG)]);
    title('BG Subtracted + Changepoint');
    legend('PL-BG','CP',thresh);
    subplot(2,1,2);
    plot(linspace(0,length(OnOff)/fps,length(OnOff)),OnOff,'k');
    ylim([0 1.5]); xlabel('Time (s)');
    title('Processed On-Off Based');
    exportgraphics(figure(4),name);
    %clear figures 
    clf(figure(3)); clf(figure(4));
    results(j,:) = [fluct,c_on.p1,gof_on.rmse,c_off.p1,gof_off.rmse];
    fprintf(fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n',results(j,:));
end
fprintf(fid,'Average values\n');
fprintf(fid,'%.2f,%.2f,%.2f,%.2f,%.2f\n',nanmean(results));
fclose(fid);
type(fileName);

%% Fit all data
all_on = [s_on clust_on]; all_off = [s_off clust_off];
[c,gof] = fit(log10(all_on(1,:)'),log10(all_on(2,:)'),'poly1');
m = sprintf('m_{on} = %.2f (%.2f)',[c.p1,gof.rmse]);
figure(5); fig = tiledlayout(1,2);
nexttile; loglog(all_on(1,:),10.^polyval([c.p1,c.p2],log10(all_on(1,:))),'b');
hold on
loglog(all_on(1,:),all_on(2,:),'ok'); 
legend(m);title('On TIme'); xlabel('log(time) (s)'); ylabel('log(PD) (1/s)');
nexttile; 
[c,gof] = fit(log10(all_off(1,:)'),log10(all_off(2,:)'),'poly1');
m =  sprintf('m_{off} = %.2f (%.2f)',[c.p1,gof.rmse]);
loglog(all_off(1,:),10.^polyval([c.p1,c.p2],log10(all_off(1,:))),'b');
hold on
loglog(all_off(1,:),all_off(2,:),'ok');
legend(m); title('Off TIme'); xlabel('log(time) (s)');
img = "Fit of all "+source+".png"; name = fullfile(F_save,img);
exportgraphics(fig,name);

%% Fit all SINGLE QD data only
clf(figure(5)); 
[c,gof] = fit(log10(s_on(1,:)'),log10(s_on(2,:)'),'poly1');
m = sprintf('m_{on} = %.2f (%.2f)',[c.p1,gof.rmse]);
figure(5); fig = tiledlayout(1,2);
nexttile; loglog(s_on(1,:),10.^polyval([c.p1,c.p2],log10(s_on(1,:))),'b');
hold on
loglog(s_on(1,:),s_on(2,:),'ok'); 
legend(m);title('On TIme'); xlabel('log(time) (s)'); ylabel('log(PD) (1/s)');
nexttile; 
[c,gof] = fit(log10(s_off(1,:)'),log10(s_off(2,:)'),'poly1');
m =  sprintf('m_{off} = %.2f (%.2f)',[c.p1,gof.rmse]);
loglog(s_off(1,:),10.^polyval([c.p1,c.p2],log10(s_off(1,:))),'b');
hold on
loglog(s_off(1,:),s_off(2,:),'ok');
legend(m); title('Off TIme'); xlabel('log(time) (s)');
img = "Fit of all Single "+source;
sgtitle(img);
name = fullfile(F_save,img+".png");
exportgraphics(fig,name);

%% Fit all CLUSTER data only
total_clusters = nansum(cluster(:,1));
if(total_clusters > 0)
    clf(figure(5)); 
    [c,gof] = fit(log10(clust_on(1,:)'),log10(clust_on(2,:)'),'poly1');
    m = sprintf('m_{on} = %.2f (%.2f)',[c.p1,gof.rmse]);
    figure(5); fig = tiledlayout(1,2);
    nexttile; loglog(clust_on(1,:),10.^polyval([c.p1,c.p2],log10(clust_on(1,:))),'b');
    hold on
    loglog(clust_on(1,:),clust_on(2,:),'ok'); 
    legend(m);title('On TIme'); xlabel('log(time) (s)'); ylabel('log(PD) (1/s)');
    nexttile; 
    [c,gof] = fit(log10(clust_off(1,:)'),log10(clust_off(2,:)'),'poly1');
    m =  sprintf('m_{off} = %.2f (%.2f)',[c.p1,gof.rmse]);
    loglog(clust_off(1,:),10.^polyval([c.p1,c.p2],log10(clust_off(1,:))),'b');
    hold on
    loglog(clust_off(1,:),clust_off(2,:),'ok');
    legend(m); title('Off TIme'); xlabel('log(time) (s)');
    img = "Fit of all Cluster "+source;
    sgtitle(img);
    name = fullfile(F_save,img+".png");
    exportgraphics(fig,name);
end

fclose('all'); 

    