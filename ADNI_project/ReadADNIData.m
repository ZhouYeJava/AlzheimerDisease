function [data_all, data_labels, RIDs, ages_all, dx_all] =...
                                        ReadADNIData(filenames,biomarkers)
% New version (ReadADNIData.m) does not rely on matlabs importdata
% function. It also relies on the data in the csv file to be organized as
% by Zhou, as listed in ReadMe.pdf (The first 6 columns must be PHASE,
% RID, VISCODE, AGE, EXAMDATE, DOB). The input 'biomarkers' is a vector of
% indices into the biomarkers listed in ReadMe.pdf (i.e. biomarkers = [7 14
% 15] will read CURRENTDIAG, RelativeHippo, MMSE)
%
% Note: Will not read Phase 2 (ADNI2, VISCODE = vXX) data! It is easily
% modified to do so  by adding cases 

%% Get data from file
nheaderlines = 1;
[data_scores, header] = getTextData(filenames,nheaderlines);

%% Make sure we have RID, VISCODE, AGE columns
if strcmp(header{2},'RID')
    RIDs = data_scores{1}(:,2);
else
    error('Incorrect file format: RID must be in the second column')
end

if strcmp(header{3},'VISCODE')
    VISCODEs = data_scores{2};
else
    error('Incorrect file format: VISCODE must be in the third column')
end

if strcmp(header{4},'AGE')
    ages = data_scores{3};
else
    error('Incorrect file format: AGE must be in the fourth column')
end

data = data_scores{5};

clear data_scores

%% Extract the biomarkers that we want
% first 6 values columns of the data are read separately
biomarkers = biomarkers - 6; 

% If we have a negative index, it should be age
inds = find(biomarkers < 0);

% Flag if age is included as a biomarker
ageflag = false;
if ~isempty(inds)
    if biomarkers(inds) ~= -2
        error(['Invalid biomarker ' header{biomarkers(inds)+6}])
    else
        ageflag = true;
    end
end
biomarkers(inds) = [];

header = header(7:end);

% Make sure we have current diagnosis
[~, bl] = ismember('CURRENTDIAG',header);
if bl == 0
    error('Missing CURRENTDIAG!')
else
    dx = data(:,bl);
end

% Get the rest of the data
ind = find(biomarkers > size(data,2),1);
if ~isempty(ind)
    error(['Biomarker ' biomarkers(ind) ' not in the data!'])
end
data = data(:,biomarkers);
header = header(biomarkers);

% Organize data by RID and visit
ages_all = nan(max(RIDs),50);
dx_all = nan(max(RIDs),50);
data_all = nan(max(RIDs),size(data,2),50);
max_ind = 0;
for i = 1:size(data,1)
    di = data(i,:);
    di(di==-6) = nan;
    rid = RIDs(i,:);
    vc = VISCODEs{i,:};
    switch vc
        case {'bl','sc'}
            ind = 1;
        case {'m03'}
            ind = 2;
        case {'m06','m6'}
            ind = 3;
        case 'm12'
            ind = 4;
        case 'm18' %mci only
            ind = 5;
        case 'm24'
            ind = 6;
        case 'm30'
            ind = 7;
        case 'm36' %n,mci only
            ind = 8;
        case 'm42'
            ind = 9;
        case 'm48' %n,mci only
            ind = 10;
        case 'm54'
            ind = 11;
        case 'm60'
            ind = 12;
        case 'm66'
            ind = 13;
        case 'm72'
            ind = 14;
        case 'm78'
            ind = 15;
        otherwise
            % Skip ADNI2
            continue
%             error(['VISCODE ' vc ' not allowed'])
    end
    if ind > max_ind
        max_ind = ind;
    end
    data_all(rid,:,ind) = di;
    ages_all(rid,ind) = ages(i);
    dx_all(rid,ind) = dx(i);
end

% Cut empty visits
data_all(:,:,(max_ind+1):end) = [];
ages_all(:,(max_ind+1):end) = [];
dx_all(:,(max_ind+1):end) = [];

% Add age biomarker 
if ageflag
    data_all = cat(2,...
        reshape(ages_all,[size(ages_all,1) 1 size(ages_all,2)]),...
        data_all);
    header = cat(2,'AGE',header);
end

RIDs = 1:max(RIDs);

data_labels = header;

    
%-------------------------------------------------------------------------%
function [data, header] = getTextData(filename,hlines)

try
    fid = fopen(filename);
catch
    error([filename ' does not exist'])
end

headerLine = fgetl(fid);
cellLine = split(headerLine, ',');
header = cellLine;

numHeaderCols = size(cellLine,2);

% Assume ordering of data is PHASE, RID, VISCODE, AGE, EXAMDATE, DOB,
%   followed by biomarker data
formatString = ['%f%f%q%f%q%q' repmat('%f', 1, numHeaderCols-6)];

data = textscan(fid,formatString,'delimiter',',',...
            'headerlines',hlines, 'CollectOutput', true);

fclose(fid);

%-------------------------------------------------------------------------%
function [cellOut, indOut] = split(fileString, delim)

cellOut = textscan(fileString,'%q','delimiter',delim,...
    'whitespace','');
indOut = strfind(fileString,delim);

cellOut = (cellOut{1})';