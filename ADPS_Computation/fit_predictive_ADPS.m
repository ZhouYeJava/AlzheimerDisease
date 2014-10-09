% Main script to fit the subject data to an already made sigmoid ADPS model
clear all
% close all

%% - Load model and data
load results/results-sigmoid-12biomarkers-16-Apr-2013

fitOpts.pos_slope=0;

opts = fitOpts;
model = fitOpts.model;
m_params = results.m_params;
s_params = results.s_params;
cov_params = results.cov_params;
data_labels = fitOpts.data_labels;

% - Last visit to use for calculating predictive ADPS
max_visit = 'bl';

%% - Get normalization factors so that the ADPS and sigmoids are on the same
%    scale as when fit using all of the data

% Normalize ADPS so that the median of normal subjects is 0 and the median
%   absolute deviation is 1
adps = repmat(s_params(:,1),[1 size(ages,2)])+repmat(s_params(:,2),[1 size(ages,2)]).*ages;
adps_mean = nanmean(adps(:)); adps_std = nanstd(adps(:));
new_adps_std = adps_std/(1.4826*mad(adps(dx==1),1));
new_adps_mean = (adps_mean-nanmedian(adps(dx==1)))/(1.4826*mad(adps(dx==1),1));

%% - Remove visits after max_visit
viscodes = {'bl','m03','m06','m12','m18','m24','m30','m36','m42','m48',...
            'm54','m60','m66','m72','m78'};
ind = find(strcmp(viscodes,max_visit));
viscodes = viscodes(1:ind);

data = data(:,:,1:ind);
ages = ages(:,1:ind);
dx_r = dx(:,1:ind);

%% - Fit to sigmoid model
[s_params,err] = fitSubjectsToSigmoid(data,ages,m_params,results.var_k,fitOpts,s_params);
err

% Normalize
s_params(:,2) = new_adps_std*s_params(:,2)/adps_std;
s_params(:,1) = new_adps_std*s_params(:,1)/adps_std + new_adps_mean - new_adps_std*adps_mean/adps_std;
adps = repmat(s_params(:,1),[1 size(ages,2)])+repmat(s_params(:,2),[1 size(ages,2)]).*ages;
m_params(:,2) = m_params(:,2)*adps_std/new_adps_std;
m_params(:,3) = new_adps_std/adps_std*(m_params(:,3)+adps_mean-new_adps_mean*adps_std/new_adps_std);

%% - Save to spreadsheet
answ = questdlg('Do you want to save the data to a spreadsheet?',...
    'Save data?','Yes','No','No');

if strcmp(answ,'Yes')
    if length(viscodes) ~= size(adps,2)
        error('number of visits for ADPS not equal to number of viscodes')
    end
    
    fn = sprintf('ADPS_2013_%s.csv',max_visit);
    
    fid = fopen(fn,'w');
    fwrite(fid,'RID');
    for i = 1:length(viscodes)
        fprintf(fid,',%s',viscodes{i});
    end
    fprintf('\r\n')
    fclose(fid);
    dlmwrite(fn,[RIDs' adps],'-append','roffset',1);
end

%% --- Unstandardize results
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

%% - Plot results
[numSubjects, numBiomarkers, numVisits] = size(data);

if numBiomarkers <= 6
    ppr = 3; %plots per row
elseif numBiomarkers < 13
    ppr = 4;
else
    ppr = 6;
end
numRows = ceil(numBiomarkers/ppr);

adps_plot = linspace(-15,15,1000);

if strcmp(model,'sigmoid')
    maxslope = -m_params(:,3);
end

yva = data; dxa = dx_r; adpsa = adps;
rp = randperm(size(data,1)); 
dxa = dxa(rp,:); adpsa = adpsa(rp,:);

figure 
for i = 1:numBiomarkers
    subplot(numRows,ppr,i)
    yvai = squeeze(data(:,i,:));
    yvai = yvai(rp,:); 
    scatter(adpsa(:),yvai(:),5,dxa(:));
    hold on
    plot(adps_plot,feval(@logisticfun,m_params(i,:),adps_plot),'Color','b','LineWidth',2)
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