function [s_params,err] = fitSubjectsToSigmoid(data,ages,m_params,var_k,opts,init_params,cov_params)

if nargin < 7
    cov_params = [];
end

[numSubjects, numBiomarkers, numVisits] = size(data);
modelfun = getModelFun(opts.model);

if opts.pos_slope
    options = optimset('Display','off','Algorithm','trust-region-reflective');
    lb = [-inf 0];
else
    options = optimset('Display','off','Algorithm','levenberg-marquardt');
    lb = [];
end

options2 = optimset('Display','off','LargeScale','off');

% Initialize
s_params = init_params;

data_i = data;
for i = 1:numBiomarkers
    di = data_i(:,i,:);
    a = m_params(i,1);
    d = m_params(i,4);
    di(di<=d) = d+.01;
    di(di>=d+a) = d+a-.01;
    data_i(:,i,:) = di;
end

wts = 1./repmat(sqrt(var_k),[1 numVisits]);

gam_X = zeros(numSubjects,numBiomarkers);
if ~isempty(cov_params)
    gam_X = covariates*cov_params';
end

err = 0;
for i = 1:numSubjects
% parfor i = 1:numSubjects
    if mod(i,10) == 0
        fprintf('%d ',i)
    end
    if mod(i,200) == 0
        fprintf('\n')
    end

    m_pi = m_params;
    m_pi(:,3) = m_pi(:,3) + gam_X(i,:)';
    
    % Get subject data
    dvi = squeeze(data(i,:,:));
    agesi = ages(i,:);
    
    if size(dvi,1) == 1
        dvi = dvi';
    end

    % fminsearch is sometimes more accurate than lsqcurvefit?
    [s_params(i,:), erri,ef] = fminsearch(@(s_p) fit_subjects(s_p,...
        dvi,agesi,m_params,modelfun,wts),s_params(i,:),options2);

%     % Weight subject data
%     dvi = dvi.*wts; nv = isnan(dvi); dvi(nv) = [];
%     
%     % Fit to model
%     [s_params(i,:), erri, ~, ef] = lsqcurvefit(@(abi,xdata) patientfun(abi,xdata,m_pi,modelfun,nv,wts),...
%                 s_params(i,:),agesi,dvi',lb,[],options);
    
    % Error of fit
    err = err + erri;
end
fprintf('\n')

function fitval = fit_subjects(sp,yvals,ages,mp,mfun,wtm)

zvals = sp(1) + ages*sp(2);
yp = feval(mfun,mp,zvals);
ferr = wtm.*(yvals - yp);
fv = isfinite(ferr);
fitval = ferr(fv)'*ferr(fv);

% --- Function to get model ---
function modelfun = getModelFun(model)

switch model
    case 'linear'
        modelfun = @linearfun;
    case 'sigmoid'
        modelfun = @logisticfun;   %@sigmoidfun
    case 'richards'
        modelfun = @richardsfun;
    case 'quadratic'
        modelfun = @quadraticfun;
    otherwise
        error('Invalid model; must be either linear or sigmoid')
end

% --- Function to evaluate fit of subject to model ---
function yvals = patientfun(a_coefs,ages,m_params,modelfun,fv,wts)
% Note jacobian would be B(y-D)*(1-(y-D)/A) and t*B(y-D)*(1-(y-D)/A)

adps = a_coefs(1) + ages*a_coefs(2);
yvals = feval(modelfun,m_params,adps);
yvals = yvals.*wts;
yvals(fv) = [];
yvals = yvals';
