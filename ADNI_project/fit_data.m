% Main script to fit the data to the ADPS model
clear all
close all

% Code can be slightly parallelized using the Parallel Computing Toolbox,
%   uncomment below to do parallel processing, and also change line 207 in
%   calculateADPS.m to parfor
s = matlabpool('size');
if s == 0
    mpl = questdlg('Matlabpool not opened. Would you like to enable it?');
    if strcmp(mpl,'Yes')
        matlabpool
    end
end

%% --- Parameters to Set ---

%--- Data file
datapath = '.\Data\';
datafilename = 'ADNI1_ADNIGO_ALL_DATA.csv';
filenames=[datapath datafilename];

%--- Select biomarkers to use - see ReadMe.pdf for indices 
biomarkers = [4 7 8 14 15 26 31 33 34 39 59 60];
% biomarkers = [14 8 15 33 34 16 26];
nonbiomarkers = []; % You can read other data for post-analysis if you want it

%--- Preprocessing parameters

% Do a regression of the data on covariates (ie regress on age, gender, etc)
ppOpts.covariates = []; % based on biomarker index in ReadMe.m
% ppOpts.covariates = [4 7 13];

% Minimum number of visits required per subject
ppOpts.minVisits = 2;

% Which biomarker is required? I often require each subject to have one
%   measurement of hippocampal volume
ppOpts.requiredBiomarker = 14; % Using ReadMe.pdf biomarker indices
% ppOpts.requiredBiomarker = 7;

% How do we handle missing data?
ppOpts.missingDataType = 1; % 1 - no change, 2 - must have all measurements, 3 - impute measurements (not implemented)

% Standardize the data to have zero mean and standard deviation one?
ppOpts.standardizeData = true;

% Exclude ABeta negative subjects
ppOpts.ABetaPos = false;

%--- Fitting algorithm parameters
% What model to fit? linear, sigmoid, sigmoidMaxMin, sigmoidMaxMinDX or sigmoid2D
% Note that we have only explored the linear and sigmoid models
%   - sigmoidMaxMin fits a sigmoid with a fixed max and min value
%     (asymptotes) (see get_biomarker_bounds.m to set these values)
%   - sigmoidMaxMinDx fixes the max and min of the sigmoid based on the
%     average of the normal and AD populations
%   - sigmoid2D is Brunos experimental 2D model
fitOpts.model = 'sigmoid';

% Weight biomarker values by variance of model fit
fitOpts.useVar = true;

% Use an autoregressive (AR(1)) model - NOT IMPLEMENTED
fitOpts.ar1model = false;

% Maximum number of fitting iterations - note there is also a convergence
%   criteria hard coded
fitOpts.MaxIter = 30; 

% Plot data at each iteration
fitOpts.doPlots = true;

% Fit data using zero mean ages (might help with stability of algorithm,
%   but I am not sure)
fitOpts.ageOffset = true;

% Fit ADPS with a correction for a covariate (s = alpha*t + beta + cov)
% fitOpts.covariates = 42;
fitOpts.covariates = [];


%% --- Read and process data ---

%--- Read data
[data, data_labels, RIDs, ages, dx] = ReadADNIData(filenames,biomarkers);

% % Fix error in trabscor data
ind = strcmp(data_labels,'TRABSCOR');
if ~isempty(ind)
    dvi = data(:,ind,:);
    dvi(dvi >= 300) = nan; % maximum time is 300
    dvi(dvi == 0) = nan; % impossible value
    data(:,ind,:) = dvi;
end

% Pre-processing covariate data
if ~isempty(ppOpts.covariates)
    cov_data = ReadADNIData(filenames,ppOpts.covariates);
    ppOpts.covariate_data = cov_data;
end

% Algorithm covariate data
if ~isempty(fitOpts.covariates)
    cov_data = ReadADNIData(filenames,fitOpts.covariates);
    fitOpts.covariate_data = cov_data;
end

% Other biomarkers
if ~isempty(nonbiomarkers)
    [nbm_data, nbm_labels] = ReadADNIData(filenames,nonbiomarkers);
else
    nbm_data = []; nbm_labels = [];
end

fitOpts.data_labels = data_labels;

%--- Preprocess data
if ~isempty(ppOpts.requiredBiomarker)
    ind = find(biomarkers == ppOpts.requiredBiomarker);
    if isempty(ind)
        error('invalid value for ppOpts.requiredBiomarker')
    end
    ppOpts.requiredBiomarker = ind;
end

[data, ages, dx, data_stats] = PreprocessADNIData(...
                    data,ages,dx,ppOpts);
                
% Trim data
inds = any(isfinite(ages),2);
data = data(inds,:,:);
ages = ages(inds,:);
dx = dx(inds,:);
RIDs = RIDs(inds);
if ~isempty(fitOpts.covariates)
    fitOpts.covariate_data = cov_data(inds,:,:);
end
if ~isempty(nonbiomarkers)
    nbm_data = nbm_data(inds,:,:);
end

% % Need a better way to set these bounds...
% if strcmp(fitOpts.model,'sigmoidMaxMin')
%     fitOpts.bounds = get_biomarker_bounds(data_labels);
%     if ppOpts.standardizeData
%         fitOpts.bounds(:,1) = (fitOpts.bounds(:,1)-data_stats(:,1))./data_stats(:,2);
%         fitOpts.bounds(:,2) = (fitOpts.bounds(:,2)-data_stats(:,1))./data_stats(:,2);
%     end
%     fitOpts.model = 'sigmoid';
% elseif strcmp(fitOpts.model,'sigmoidMaxMinDX')
%     dx_all = squeeze(nbm_data(:,2,:));
%     fitOpts.bounds = get_biomarker_bounds_DX(data,dx_all);
%     fitOpts.model = 'sigmoid';
% end


%% --- Fit model to data ---

%--- Initialize using linear fit
fitOptsLinear.useVar = false;
fitOptsLinear.ar1model = false;
fitOptsLinear.maxIter = 100;
fitOptsLinear.doPlots = true; 
fitOptsLinear.data_labels = data_labels;
fitOptsLinear.model = 'linear';
fitOptsLinear.s_init = [zeros(size(data,1),1) ones(size(data,1),1)];
fitOptsLinear.ageOffset = true;
[mp, s_params var_k cov_params conv_vals] = calculateADPS(data,ages,fitOptsLinear);

csvwrite('s_params_after_linear_fit.csv',s_params)

%--- Fit model to data
fitOpts.s_init = s_params;
fitOpts.pos_slope = false;
if ~strcmp(fitOpts.model,'sigmoid2d')
    [m_params s_params var_k cov_params conv_vals] = calculateADPS(data,ages,fitOpts);
else
    [m_params s_params var_k cov_params] = calculateADPS2d(data,ages,fitOpts);
end

%--- Make sure hippocampus has negative slope
if strcmp(fitOpts.model,'sigmoid')
    m_params = convertLogisticParams(m_params);
    % NOTE: I assume hippocampus is the first biomarker here!
    ind = strcmp(data_labels,'RelativeHippo');
    if isempty(ind)
        warning('No hippocampus value to ensure the ADPS is going in the correct direction (N->AD)')
    else
        if m_params(ind,2) > 0
            m_params(:,2:3) = -m_params(:,2:3);
            s_params = -s_params;
        end
    end
end

% %--- Unstandardize results
% if ~isempty(data_stats)
%     dvm = data_stats(:,1); dvstd = data_stats(:,2);
%     % unstandardize model parameters
%     switch fitOpts.model
%         case 'linear'
%             m_params(:,2) = m_params(:,2).*dvstd;
%             m_params(:,1) = m_params(:,1).*dvstd + dvm;
%         case 'sigmoid'
%             m_params(:,1) = m_params(:,1).*dvstd;
%             m_params(:,4) = m_params(:,4).*dvstd + dvm;
%         case 'sigmoid2d'
%             m_params(:,1) = m_params(:,1).*dvstd;
%             m_params(:,6) = m_params(:,6).*dvstd + dvm;
%     end
%     % unstandardize data
%     data = StandardizeData(data,dvm,dvstd); 
% end

results.m_params = m_params;
results.s_params = s_params;
results.var_k = var_k;
results.cov_params = cov_params;

% adps = repmat(s_params(:,1),[1 size(ages,2)])+repmat(s_params(:,2),[1 size(ages,2)]).*ages;
% [adps, s_params, m_params, cov_params] = normalizeADPS(adps,dx,s_params,m_params,cov_params,fitOpts.model,2);


%% --- Save data and results --
if isfield(fitOpts,'s_init')
    fitOpts = rmfield(fitOpts,'s_init');
end
btn = questdlg('Save results?');
if strcmp(btn,'Yes')
    if ~exist('results','dir')
        mkdir('results')
    end
    savefile = ['results/results-' fitOpts.model '-' num2str(numel(biomarkers)) 'biomarkers-' date];
    save(savefile,'data','ages','dx','RIDs','results','biomarkers','filenames',...
        'ppOpts','fitOpts','nbm_data','nbm_labels','data_stats');
end

