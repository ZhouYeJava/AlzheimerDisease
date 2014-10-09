function [data, ages, dx, data_stats] = PreprocessADNIData(...
                        data,ages,dx,opts)
% Preprocessing steps (as dictated by opts):
%   1) Adjust the data for covariates (using regression)
%       - Covariates are those in the opts.covariate_data
%   2) Remove visits without diagnosis
%   3) Remove visits without a subject age
%   4) Remove visits without required biomarker
%   5) Modify the data to account for missing data, there are multiple
%       options for accounting for missing data:
%           1 - leave as is
%           2 - complete data case
%           3 - imputation (to be implemented? regression?)
%           4 - must have x% of data (to be removed?)
%   6) Remove subjects without the minimum number of non-empty visits
%   7) Standardize the data to mean 0 and std 1 (z-scores)
%   8) Remove age and diagnosis without any data 
%       (ie remove nonbiomarker data when there is no actual data there,
%       for instance, there could be an age listed without any data)

if nargin < 3 % technically don't need nbm data?
    error('Invalid number of input arguments');
elseif nargin == 3
    opts = struct([]);
end
opts = getopts(opts); % fill in defaults

%% --- Adjust for covariates ---
% note that visits without all covariates are implicitly removed by the
%   regression
% also note, we only do the regression on normal subjects, but apply the
%   adjustment to all subjects
covariates = opts.covariates;
if ~isempty(covariates)
    covariates = opts.covariate_data;
    if isempty(covariates)
        error('Covariate data missing for preprocessing step!')
    end
    
    % Do regression
    [data adjParams] = adjustForCovariates(data,covariates,dx);
end

%% --- Remove visits without a diagnosis given ---
for i = 1:size(data,1)
    dvi = squeeze(data(i,:,:));
    dxi = dx(i,:);
    dvi(:,isnan(dxi)) = nan;
    data(i,:,:) = dvi;
end


%% --- Remove visits without an age given ---
for i = 1:size(data,1)
    dvi = squeeze(data(i,:,:));
    agei = ages(i,:);
    dvi(:,isnan(agei) | agei < 10) = nan;
    data(i,:,:) = dvi;
end    


%% --- Remove data without required biomarker ---
if ~isempty(opts.requiredBiomarker)
    rmb = opts.requiredBiomarker;
    for i = 1:size(data,1)
        dvi = squeeze(data(i,:,:));
        nv = isnan(dvi);
        
        % remove *visits* without required biomarker
        dvi(:,any(nv(rmb,:),1)) = nan;
        data(i,:,:) = dvi;
        
%         % remove *subjects* without required biomarker
%         if all(nv(rmb,:)==1)
%             data(i,:,:) = nan;
%         end
    end
end

%% --- Remove ABeta negative subjects ---
if ~(isempty(opts.ABetaPos) || opts.ABetaPos == false)
    % Get covariate data from nbm_data
    [~, cl] = ismember('ABETA',nbm_labels);
    if cl==0
        error('Missing ABeta data from nonbiomarker data')
    end
    abeta = squeeze(nbm_data(:,cl,:));
    
    abeta_pos = abeta < 192;
    abeta_neg = abeta >= 192;
    
    anyabeta = any(abeta_pos,2);
    data(~anyabeta,:,:) = nan;
end


%% --- Adjust for missing data ---
% 1 = no change
% 2 = must have all data
% 3 = imputation (to be implemented?)
if opts.missingDataType ~= 1
    missingDataType = opts.missingDataType;
    for i = 1:size(data,1)
        dvi = squeeze(data(i,:,:));

        if missingDataType == 2
            % must have all data at each visit, ie. remove visits without
            %   all biomarkers
            nv = isnan(dvi);
            dvi(:,any(nv,1)) = nan;
            data(i,:,:) = dvi;
        elseif missingDataType == 3
            % imputation
            
        end
    end
end


%% --- Remove subjects without minimum number of visits ---
minVisits = opts.minVisits;

if minVisits > 1
    for i = 1:size(data,1)
        dvi = squeeze(data(i,:,:));
        
        nv = sum(any(isfinite(dvi),1));
        
        if nv < minVisits && nv > 0
            data(i,:,:) = nan;
        end
    end
end


%% --- Standardize the data to mean 0, std 1 ---
% Should this be done before data is removed (steps 1-6)?
if opts.standardizeData    
    [data, data_mean, data_std] = StandardizeData(data);
    data_stats = [data_mean data_std];
else
    data_stats = [];
end


%% --- Remove visits (age and dx) without data --
% Force age and diagnosis to be NaN if there is no data for that visit
for i = 1:size(data,1)
    dvi = squeeze(data(i,:,:));
    allgone = all(isnan(dvi),1);
    ages(i,allgone) = nan;
    dx(i,allgone) = nan;
end



function newopts = getopts(opts)
% Function to get default options

% Defaults
options = struct('covariates', [], 'minVisits', 2,...
    'requiredBiomarker', [], 'missingDataType', 1,...
    'standardizeData', true, 'ABetaPos', false, 'covariate_data', []);

names = fieldnames(options);
numNames = numel(names);

argNames = fieldnames(opts);
newopts = options;

for i = 1:numNames
    name = names{i};
    if any(strcmp(name,argNames))
        val = opts.(name);
        if ~isempty(val)
%             if ischar(val)
%                 val = lower(deblank(val));
%             end
%             checkparam(name,val);
            newopts.(name) = val;
        end
    end
end

function [data, adjParams] = adjustForCovariates(data,covars,dx)
% Adjust measurements for change in covariates (age,education,etc...)
%   covars should be an nxmxc, n is # subjects, m is # visits, c is #
%   covariates
% If there is a diagnosis, we do a regression only on the normal subjects
%   and apply result to all subjects, otherwise, do for all subjects

% Assume minimum diagnosis value is control
if isempty(dx)
	dx = ones(size(data,1),size(data,3));
end
dx_min = min(dx(:));
dx_control = dx==dx_min;
% dx_control = true(size(dx)); % uncomment to use all subjects in regression

if size(covars,2) ~= size(dx,2)
    covars = permute(covars,[1 3 2]);
end

covars = reshape(covars,[],size(covars,3));
covarsc = covars(dx_control(:),:);

[numSubjects, numBiomarkers, numVisits] = size(data);

adjParams = zeros(numBiomarkers,size(covarsc,2)+1);
for i = 1:numBiomarkers
    data_i = squeeze(data(:,i,:));
    
    d_c = data_i(dx_control);  

    % Do regression on all data vs covariates
    [adjParams(i,:),bint,r,rint,stats] = regress(d_c,[ones(size(d_c)) covarsc]);

    % Adjusted data is residuals plus mean of data
    fv = isfinite(r);
    d_m = mean(d_c(fv));
    d_m = 0;
    data_i = data_i(:) - [ones(size(covars,1),1) covars]*adjParams(i,:)' + d_m;
    data_i = reshape(data_i,numSubjects,numVisits);
    data(:,i,:) = data_i;
end
