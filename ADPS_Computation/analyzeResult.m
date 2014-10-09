clear all
close all


%% --- Load Data ---

load results/results-sigmoid-12biomarkers-18-Apr-2013
opts = fitOpts;
model = fitOpts.model;
m_params = results.m_params;
s_params = results.s_params;
cov_params = results.cov_params;
data_labels = fitOpts.data_labels;

switch model
    case 'linear'
        modelfun = @linearfun;
    case 'sigmoid'
        modelfun = @logisticfun;   %@sigmoidfun
    otherwise
        error('Invalid model; must be either linear or sigmoid')
end

[numSubjects numBiomarkers numVisits] = size(data);

%% --- Unstandardize results --
if ~isempty(data_stats)
    dvm = data_stats(:,1); dvstd = data_stats(:,2);
    % unstandardize model parameters
    switch fitOpts.model
        case 'linear'
            m_params(:,2) = m_params(:,2).*dvstd;
            m_params(:,1) = m_params(:,1).*dvstd + dvm;
        case 'sigmoid'
            m_params(:,1) = m_params(:,1).*dvstd;
            m_params(:,4) = m_params(:,4).*dvstd + dvm;
        case 'sigmoid2d'
            m_params(:,1) = m_params(:,1).*dvstd;
            m_params(:,6) = m_params(:,6).*dvstd + dvm;
    end
    % unstandardize data
    data = StandardizeData(data,dvm,dvstd); 
end

%% --- Normalize data ---

% Calculate ADPS
adps = repmat(s_params(:,1),[1 size(ages,2)])+repmat(s_params(:,2),[1 size(ages,2)]).*ages;
% [adps, s_params, m_params, cov_params] = normalizeADPS(adps,dx,s_params,m_params,cov_params,model,2);

% Normalize ADPS so that the median of normal subjects is 0 and the median
%   absolute deviation is 1 
adps_mean = nanmean(adps(:)); adps_std = nanstd(adps(:));
new_adps_std = adps_std/(1.4826*mad(adps(dx==1),1));
new_adps_mean = (adps_mean-nanmedian(adps(dx==1)))/(1.4826*mad(adps(dx==1),1));
% Normalize ADPS so that the mean of normal subjects is 0 and the standard
%   deviation is 1
% new_adps_std = adps_std/nanstd(adps(dx==1));
% new_adps_mean = (adps_mean-nanmean(adps(dx==1)))/nanstd(adps(dx==1));
adps = new_adps_std*(adps - adps_mean)/adps_std+new_adps_mean;

% Normalize the subject and model parameters
s_params(:,2) = new_adps_std*s_params(:,2)/adps_std;
s_params(:,1) = new_adps_std*s_params(:,1)/adps_std + new_adps_mean - new_adps_std*adps_mean/adps_std;
switch model
    case 'linear'
        m_params(:,2) = m_params(:,2)*adps_std/new_adps_std;
        m_params(:,1) = m_params(:,1)+m_params(:,2)*(new_adps_std*adps_mean/adps_std-new_adps_mean);
    case 'sigmoid'
        m_params(:,2) = m_params(:,2)*adps_std/new_adps_std;
        m_params(:,3) = new_adps_std/adps_std*(m_params(:,3)+adps_mean-new_adps_mean*adps_std/new_adps_std);
end


% ad_agesp = linspace(min(adps(:))-2,max(adps(:))+2,100000);

% The range of ADPS we want to plot over (dont visualize outliers)
adps_m = nanmean(adps,2);
qv = quantile(adps_m,[0.025 0.975]);
% Add a buffer to the inter-quantile range
qvr = 1*(qv(2)-qv(1));        
qv(1) = qv(1) - qvr; qv(2) = qv(2) + qvr;
adps_plot = linspace(qv(1),qv(2),10000);

ve = explainedVariance(data,m_params,adps,[]);

%% --- Save to spreadsheet? ---

viscodes = {'bl','m03','m06','m12','m18','m24','m30','m36','m42','m48',...
            'm54','m60','m66','m72','m78'};
ans = questdlg('Do you want to save the data to a spreadsheet?',...
    'Save data?','Yes','No','No');

if strcmp(ans,'Yes')
    if length(viscodes) ~= size(adps,2)
        error('number of visits for ADPS not equal to number of viscodes')
    end
    
    fid = fopen('ADPS_2013.csv','w');
    fwrite(fid,'RID');
    for i = 1:length(viscodes)
        fprintf(fid,',%s',viscodes{i});
    end
    fprintf('\r\n')
    fclose(fid)
    dlmwrite('ADPS_2013.csv',[RIDs' adps],'-append','roffset',1);
end

%% --- Plot first 6 biomarkers in 3D parametric space ---

% What does our line look like in parametric space of the first 3
% biomarkers?
if numBiomarkers > 2
    figure,
    % Plot the data
    scatter3(reshape(data(:,1,:),1,[]),reshape(data(:,2,:),1,[]),reshape(data(:,3,:),1,[]),[],dx(:))
    xlabel(data_labels{1}), ylabel(data_labels{2}), zlabel(data_labels{3})
    
    % Plot the sigmoid parametric curve
    adps1 = linspace(min(adps(:)),max(adps(:)),10000);    
    pts3da = [feval(modelfun,m_params(1,:),adps1);...
              feval(modelfun,m_params(2,:),adps1);...
              feval(modelfun,m_params(3,:),adps1)];
    hold on
    plot3(pts3da(1,:),pts3da(2,:),pts3da(3,:),'k','LineWidth',2);
          
    % Plot equi-ADPS points (ie. the distance between these points
    %   represents the same increment of ADPS values)
    adps2 = linspace(min(adps(:)),max(adps(:)),20);
    pts3db = [feval(modelfun,m_params(1,:),adps2);...
              feval(modelfun,m_params(2,:),adps2);...
              feval(modelfun,m_params(3,:),adps2)];    
    plot3(pts3db(1,:),pts3db(2,:),pts3db(3,:),'b.','MarkerSize',30);
    
    colormap([0 0 0; 1 0 0; 0 1 0]);
    if numBiomarkers >= 6
    figure,
        scatter3(reshape(data(:,4,:),1,[]),reshape(data(:,5,:),1,[]),reshape(data(:,6,:),1,[]),[],dx(:))
        xlabel(data_labels{4}), ylabel(data_labels{5}), zlabel(data_labels{6})
        pts3da = [feval(modelfun,m_params(4,:),adps1);...
                  feval(modelfun,m_params(5,:),adps1);...
                  feval(modelfun,m_params(6,:),adps1)];
        pts3db = [feval(modelfun,m_params(4,:),adps2);...
                  feval(modelfun,m_params(5,:),adps2);...
                  feval(modelfun,m_params(6,:),adps2)];
        hold on
        plot3(pts3da(1,:),pts3da(2,:),pts3da(3,:),'k','LineWidth',2);
        plot3(pts3db(1,1),pts3db(2,1),pts3db(3,1),'c.','MarkerSize',30);
        plot3(pts3db(1,end),pts3db(2,end),pts3db(3,end),'m.','MarkerSize',30);
        plot3(pts3db(1,2:end-1),pts3db(2,2:end-1),pts3db(3,2:end-1),'b.','MarkerSize',30);
        colormap([0 0 0; 1 0 0; 0 1 0]);
    end
end

% What does our line look like in parametric space of the first 2
% biomarkers? Plot all pairwise biomarkers
if numBiomarkers < 8
    figure
    adps1 = linspace(min(adps(:)),max(adps(:)),10000);
    adps2 = linspace(min(adps(:)),max(adps(:)),10);
    for i = 1:numBiomarkers-1
        for j = (i+1):numBiomarkers
            subplot(numBiomarkers-1,numBiomarkers-1,sub2ind([numBiomarkers-1,numBiomarkers-1],j-1,i))
            scatter(reshape(data(:,i,:),1,[]),reshape(data(:,j,:),1,[]),[],dx(:))
            xlabel(data_labels{i}), ylabel(data_labels{j})
            pts2da = [feval(modelfun,m_params(i,:),adps1);...
                      feval(modelfun,m_params(j,:),adps1)];
            pts2db = [feval(modelfun,m_params(i,:),adps2);...
                      feval(modelfun,m_params(j,:),adps2)];
            hold on
            plot(pts2da(1,:),pts2da(2,:),'k','LineWidth',3);
            plot(pts2db(1,:),pts2db(2,:),'b.','MarkerSize',20);
        end
    end
    colormap([0 0 0; 1 0 0; 0 1 0]);
end

%% --- Plot speed vs adps colored by baseline diagnosis ---
% x-axis is average ADPS for that subject
adps_mid = zeros(numSubjects,1);
for i = 1:numSubjects
    adpsi = adps(i,:);
    adpsim = nanmean(adpsi);
    adps_mid(i) = adpsim;
end

figure
% scatter(adps(:,1),s_params(:,2),[],dx(:,1)); % baseline ADPS
nv1=(adps_mid<=qv(2))&(adps_mid>=qv(1));
scatter(adps_mid(nv1),s_params(nv1,2),[],dx(nv1,1)); % average ADPS
looseAxis;
hline(0,'k:');
xlabel('\bfADPS'), ylabel('\bfspeed')
colormap([0 0 0; 1 0 0; 0 1 0]);
% xlswrite('ad_speed',[adps(:,1),adps_mid,s_params(:,2),dx(:,1)])

%% --- Plot fitted model ---
figure
if numBiomarkers <= 6
    ppr = 3; %plots per row
elseif numBiomarkers < 13
    ppr = 4;
else
    ppr = 6;
end
numRows = ceil(numBiomarkers/ppr);

for i = 1:numBiomarkers
    subplot(numRows,ppr,i)
    yvi = squeeze(data(:,i,:));    
    plot(adps(:),yvi(:),'co','MarkerSize',2.5,'LineStyle','none')
    xlabel('\bfADPS')    

    hold on, plot(adps_plot,feval(modelfun,m_params(i,:),adps_plot),'Color','k','LineWidth',2), hold off
    
    a = axis;
    a(1:2) = [adps_plot(1),adps_plot(end)];
    axis(a);
    
    title(data_labels{i},'Interpreter','none','fontweight','b');
end

%% --- Plot each biomarker vs age color coded by diagnosis ---
yva = data; dxa = dx; 
agesa = ages;
figure
rp = randperm(size(yva,1)); 
dxa = dxa(rp,:);agesa = agesa(rp,:);
for i = 1:numBiomarkers
    subplot(numRows,ppr,i)
    yvai = squeeze(yva(:,i,:));
    yvai = yvai(rp,:); 
    scatter(agesa(:),yvai(:),5,dxa(:));
    xlabel('\bfage'),title(data_labels{i},'Interpreter','none','fontweight','b')
    set(gca,'TickDir','out')
    looseAxis;
end
colormap([0 0 0; 1 0 0; 0 1 0]);

%% --- Plot each biomarker vs ADPS color coded by diagnosis ---
if strcmp(model,'sigmoid')
    maxslope = -m_params(:,3);
end
adagesa = adps;
figure 
adagesa = adagesa(rp,:);
bds = [2000 12000; -10 190; -1.5 32; -20 400;  25 340; -2 19; -1.5 16.5];
for i = 1:numBiomarkers
    subplot(numRows,ppr,i)
    yvi = squeeze(data(:,i,:));
    yvai = squeeze(yva(:,i,:));
    yvai = yvai(rp,:); 
    scatter(adagesa(:),yvai(:),5,dxa(:));
    hold on
    plot(adps_plot,feval(modelfun,m_params(i,:),adps_plot),'Color','b','LineWidth',2)
    hold off
    looseAxis;
    a = axis;
    a(1:2) = [adps_plot(1),adps_plot(end)];
    axis(a);
    
    xlabel('\bfADPS'),title(data_labels{i},'Interpreter','none','fontweight','b')
    set(gca,'TickDir','out')
end
colormap([0 0 0; 1 0 0; 0 1 0]);
set(gcf,'Position',[184 113 1369 565])
   
%% --- Plot Cliff Jack Curves ---
figure
plotCliffJackCurve(m_params,adps_plot,data_labels);
vline(min(adps(:)),'k:'), vline(max(adps(:)),'k:')
% axis([-10.5 20.5 -0.1000 1.1000])
% hold on, plot([-10.5 20.5],[.5 .5],':','Color',[150 150 150]/255),hold off

% Lets do it again, this time try to correspond better to the paper
% hippoc(light green), adas(dark green), mmse(dark green), tau(blue),
% abeta(red), cdrsb(dark green), ravlt(purple)  
colors = distinguishable_colors(numBiomarkers);
if numBiomarkers == 7
    figure
    colors = [0.2 0.9 0;
                  0 0.5 0;
                  0 0.5 0;
                  0 0 0.9;
                  1 0 0;
                  0 0.5 0;
                  0.5 0 0.8];
    linestyles = {'-','--',':','-','-','-','-'};
    plotCliffJackCurve(m_params,adps_plot,data_labels,colors,linestyles);
    
    mci_mean = nanmean(adps(dx==2)); mci_std = nanstd(adps(dx==2));
    vline(mci_mean-mci_std,'k--'), vline(mci_mean+mci_std,'k--')
    set(gcf,'Position',[152 474 1473 504])
    a = axis; a(1:2) = [-20 20]; axis(a)
end

%% --- Plot cliff jack curve normalized in a different way ---
figure
% 1 and 99 quantile of ADPS
qle = quantile(adps(:),[0.01 0.99]);
adps_qle = linspace(qle(1),qle(2),1000);
adps_qle2 = [linspace(qle(1)-2,qle(1),100)];
adps_qle3 = [linspace(qle(2),qle(2)+2,100)];
h = zeros(numBiomarkers,1);
for i = 1:numBiomarkers
        % get sigmoid values
        mpi = m_params(i,:);
        mfe = feval(modelfun,mpi,adps_qle);
        % scale 0 to 1
        minv = min(mfe); maxv = max(mfe);
        mfe = (mfe - minv)/(maxv-minv);
        if mfe(end) < mfe(1)
            mfe = 1-mfe;
        end
        
        % get end values
        mfe2 = feval(modelfun,mpi,adps_qle2);
        mfe2 = (mfe2 - minv)/(maxv-minv);

        mfe3 = feval(modelfun,mpi,adps_qle3);
        mfe3 = (mfe3 - minv)/(maxv-minv);
        if mfe3(1) < mfe2(1)
            mfe2 = 1-mfe2;
            mfe3 = 1-mfe3;
        end
        
        hold on
        h(i) = plot(adps_qle,mfe,'LineWidth',2,'Color',colors(i,:));
        ncm = colors(i,:)+100/255; ncm(ncm>1) = 1; 
        plot(adps_qle2,mfe2,'LineWidth',2,'Color',ncm);
        plot(adps_qle3,mfe3,'LineWidth',2,'Color',ncm);
end
h = legend(h,data_labels,'Location','NorthEastOutside');
set(h,'Interpreter','none')
looseAxis; ax = axis; ax(3:4) = [-0.15 1.15]; axis(ax);
ylabel('\bfBiomarker Magnitude')
xlabel('\bfADPS')

%% --- Spaghetti plot ---
figure,
hh = zeros(6,1);
for i = 1:numBiomarkers
%     subplot(2,4,i)
    subplot(numRows,ppr,i)
%     subplot(3,3,i)
    hold all
    for j = 1:size(data,1)
        dij = squeeze(data(j,i,:));
        nv = isfinite(dij);
        dxj = dx(j,nv);
        
        if any(nv)
            if all(dxj == 1)
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(1) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','k');
            elseif all(dxj == 2)
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(2) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','r');
            elseif all(dxj == 3)
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(3) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','g');
            elseif all(dxj== 1 | dxj == 2)
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(4) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','c');
    %             cl = color_line(adps0(j,nv),dij(nv),dxj/3,'LineWidth',1.5);
            elseif all(dxj== 2 | dxj == 3)
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(5) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','b');
            else
                plot(adps(j,nv),dij(nv),'LineWidth',2,'Color','w');
                hh(6) = plot(adps(j,nv),dij(nv),'LineWidth',1.5,'Color','m');
            end
        end        
    end
    
    xlabel('\bfADPS')
    plot(adps_plot,logisticfun(m_params(i,:),adps_plot),'color','k','linewidth',2)
    plot(adps_plot,logisticfun(m_params(i,:),adps_plot),'color',[150 150 150]/255,'linewidth',1.5)
    
    title(data_labels{i},'interpreter','none','fontweight','bold')
    looseAxis;
    
    a = axis;
    a(1:2) = [adps_plot(1),adps_plot(end)];
    axis(a);
    
    set(gca,'TickDir','out')
end

% h = legend(hh,'N \rightarrow N','MCI \rightarrow MCI','AD \rightarrow AD','N \leftrightarrow MCI','MCI \leftrightarrow AD','N \leftrightarrow AD','Location','East');