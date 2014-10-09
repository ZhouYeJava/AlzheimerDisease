function [data_stdze data_mean data_stdev] = StandardizeData(data,data_mean,data_stdev)
% Standardize the data such that it is zero mean, standard deviation 1. If
% a vector of mean and stdev values are input, the inverse is applied, ie,
% try to recover original data by 'inverse standardization'

data_stdze = data;
nb = size(data,2);

if nargin < 2
    data_mean = zeros(nb,1); data_stdev = data_mean;
    for i = 1:nb
        d_i = squeeze(data(:,i,:));    
        data_mean(i) = nanmean(d_i(:));
        data_stdev(i) = nanstd(d_i(:));
        d_i = (d_i-data_mean(i))/data_stdev(i);
        data_stdze(:,i,:) = d_i;
    end
elseif nargin == 3
    for i = 1:nb
        d_i = squeeze(data(:,i,:));
        d_i = d_i*data_stdev(i)+data_mean(i);
        data_stdze(:,i,:) = d_i;
    end
else
    error('Invalid number of input arguments!')
end