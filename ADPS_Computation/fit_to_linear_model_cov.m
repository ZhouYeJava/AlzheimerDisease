function [ab e] = fit_to_linear_model_cov(yvals,Xt,covmat,u_coefs,var_k)
% Fit the data to a linear ADPS model
%
% We want to estimate alpha and beta in the model
%   y_ijk = u_1k + u_2k*(alpha_i*t_ij + beta_i) 
% where i, j, k index subject, visit and biomarker respectively
% We solve for alpha and beta by doing a regression
%   (y_ijk - u_1k)/u_2k = alpha_i*t_ij + beta_i
% Since we are dividing by u_2k, the variance of the noise is also divided
%   by u_2k ((y_ijk - u_1k)/u_2k = alpha_i*t_ij + beta_i + n_ijk/u_2k)
%
% Note the possibility of inputting a JxJ covariance matrix for the data
%   for AR models (ie. visits for the same subject are correlated)

if isempty(covmat)
    covmat = eye(size(yvals,2));
end

[numSubjects numBiomarkers numAges] = size(yvals);
datalength = numAges;

yt = zeros(datalength*numBiomarkers,numSubjects);
cm = covmat;

% Precalculate y vectors and covariance matrix
for j = 1:numBiomarkers
    y_j = (yvals(:,j,:) - u_coefs(j,1))/u_coefs(j,2);
    yt((j-1)*datalength+1:j*datalength,:) = squeeze(y_j)';    
    
    cm(j,:) = cm(j,:)/u_coefs(j,2);
    cm(:,j) = cm(:,j)/u_coefs(j,2);    
end
cm = cm.*diag(var_k);

cmk = kron(cm,eye(datalength));

% Calculate number of visits per biomarker per subject, use for weighting
nvisits = zeros(numSubjects,numBiomarkers);
for i = 1:numSubjects
    y_if = isfinite(squeeze(yvals(i,:,:)));
    nvisits(i,:) = sum(y_if,2)';
end

e = 0;
ab = zeros(numSubjects,2);
for i = 1:numSubjects

    X = Xt(:,:,i);
    y = yt(:,i);

    % Find ab by least squares fit
    yn = find(~isnan(y));
    Xyn = X(yn,:);

    if all(Xyn(:,2) == Xyn(1,2))
        % Single visit, multiple data points
        cmi = inv(cmk(yn,yn));
        yyn = y(yn);
        wts = 1./diag(cmk(yn,yn)); wts = wts/sum(wts);
        ab(i,2) = 0;
        ab(i,1) = yyn'*wts;
        
        er = (y(yn)-Xyn*ab(i,:)')'*cmi*(y(yn)-Xyn*ab(i,:)');
        e = e + er;
        
%         [abi, ~, er] = lscov(X(yn,:),yyn,cmk(yn,yn));
    else
        if length(yn) >= 2
                cmi = inv(cmk(yn,yn));

                ab(i,:) = (Xyn'*cmi*Xyn)\(Xyn'*cmi*y(yn)); 
    %             [ab(i,:), ~, er] = lscov(X(yn,:),y(yn),cmk(yn,yn));
                if ~any(isfinite(ab(i,:)))
    %                 [ab(i,:),~,er] = lscov(X(yn,:),y(yn),cmk(yn,yn));
                    ab(i,:) = NaN;
                    er = 0;
                else
    %                 er = (y(yn)-Xyn*ab(i,:)')'*cmi*(y(yn)-Xyn*ab(i,:)')/length(yn);
                    er = (y(yn)-Xyn*ab(i,:)')'*cmi*(y(yn)-Xyn*ab(i,:)');
    %                 er = (y(yn)-Xyn*ab(i,:)')'*cmi*(y(yn)-Xyn*ab(i,:)')/(length(yn)-2); %not MLE
                end
                e = e + er;
        else
            ab(i,:) = NaN;
        end
    end
end