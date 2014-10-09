function [ve vr vbm mv A] = explainedVariance(data,m_params,adps,covdata)
% Calculate the amount of explained variance of the fitted model

modelfun = @logisticfun;

% average over each biomarker
numBiomarkers = size(data,2);
A = zeros(1,numBiomarkers);
y_km = zeros(1,numBiomarkers);
for i = 1:numBiomarkers
    y_k = squeeze(data(:,i,:));
    A(i) = sum(isfinite(y_k(:)));
    y_km(i) = 1./A(i)*nansum(y_k(:));
end

% variance over each biomarker
t2 = zeros(1,numBiomarkers);
for i = 1:numBiomarkers
    y_k = squeeze(data(:,i,:));
    y_kzm2 = (y_k - y_km(i)).^2;
    t2(i) = 1./A(i)*nansum(y_kzm2(:));
end

% estimated values
y_hat = zeros(size(data));
for i = 1:size(data,1)
    for j = 1:size(data,2)
        for k = 1:size(data,3)
            if ~isempty(covdata)
                y_hat(i,j,k) = feval(modelfun,m_params(j,:),adps(i,k)+covdata(i,j));
            else
                y_hat(i,j,k) = feval(modelfun,m_params(j,:),adps(i,k));
            end
        end
    end
end

% remaining variance
v = 0;
y_sub = data-y_hat;
for i = 1:numBiomarkers
    y_sk2 = squeeze(y_sub(:,i,:)).^2;
    v = v + 1/(t2(i)*A(i))*nansum(y_sk2(:));    
end

%explained variance
ve = 1 - v/7;

vr = v;
vbm = t2;
mv = y_km;
