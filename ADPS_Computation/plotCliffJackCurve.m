function h = plotCliffJackCurve(m_params,adps,data_labels,colors,linestyles,axislabels)
% Make a plot to look similar to the infamous 'Cliff Jack Curves'
% See reference (Figure 2): 
% Jack, Jr., C. R., Knopman, D. S., Jagust, W. J., Shaw, L. M., Aisen, P.
%   S., Weiner, M. W., Petersen, R. C., Trojanowski, J. Q., 2010.
%   Hypothetical model of dynamic biomarkers of the Alzheimer's
%   pathological cascade. The Lancet Neurology 9 (1), 119-128.

numBiomarkers = size(m_params,1);

if nargin < 6
    axislabels = [];
end
if nargin < 5
    linestyles = repmat({'-'},1,numBiomarkers);
end
if nargin < 4
    colors = distinguishable_colors(numBiomarkers);
end

if isempty(colors), colors = distinguishable_colors(numBiomarkers); end
if isempty(linestyles), linestyles = repmat({'-'},1,numBiomarkers); end

% Change parameters so sigmoid goes from 0 to 1
m_params = convertLogisticParams(m_params);

hl = zeros(numBiomarkers,1);

% Plot each curve
for i = 1:numBiomarkers
    % scale data and shift appropriately
    mpi = m_params(i,:);
    mpi(4) = 0; mpi(1) = 1; if mpi(2)<0, mpi(2) = -mpi(2);end

    hold on
    plot(adps,feval(@logisticfun,mpi,adps),'LineWidth',2,'Color',colors(i,:),'LineStyle',linestyles{i})
end
if ~isempty(data_labels)
    h = legend(data_labels,'Location','NorthEastOutside');
    set(h,'Interpreter','none')
end
axis([min(adps) max(adps) -0.1 1.1])
% vline(min(adps(:)),'k:'), vline(max(adps(:)),'k:')
if isempty(axislabels)
    ylabel('\bfBiomarker Magnitude')
    xlabel('\bfADPS')
else
    ylabel(axislabels{1})
    xlabel(axislabels{2})
end
set(gca,'YTick',[0 1],'YTickLabel',{'Normal','Abnormal'},'TickDir','out');