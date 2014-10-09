function [m_params s_params var_k cov_params conv_vals] = calculateADPS(data,ages,opts)
%       - optimize z = alpha*(t-t_bar) + beta such that 
%           z = alpha*t + beta' where beta' = beta - alpha*t_bar

% If you want to add a model, need to edit 'getModelFun' and
%   'initializeModel' functions at the bottom to include it

% To Do:
%   - better plot options
%   - better handling of new models


%% --- Check inputs ---
if nargin < 2
    error('Invalid number of input arguments');
elseif nargin < 3
    opts = struct([]);
end
opts = getopts(opts); % fill in missing fields with defaults


%% --- Parameters ---
model = opts.model;
[numSubjects numBiomarkers, numVisits] = size(data);
modelfun = getModelFun(model);
bounds = opts.bounds;
ageoffset = opts.ageOffset;


%% --- Initialization stuff ---

%--- Repmat covariate data
useCovariates = false; cov_params = [];
if ~isempty(opts.covariates)
    useCovariates = true;
    covariates = opts.covariates;
%     cov_params = zeros(numBiomarkers,size(covariates,2));
    cov_params = zeros(1,size(covariates,2));
end

%--- adjust ages to zero mean for optimization
if ageoffset
    age_mean = nanmean(ages,2);
    ages_zm = bsxfun(@minus,ages,age_mean);
else
    ages_zm = ages;
end

%--- Initialization of subject parameters
if isempty(opts.s_init) && ~strcmp(model,'linear')
    % Initialize by PCA
    dva = shiftdim(data,1); dva = dva(:,:)';
    [~, ~, ~, X, ~] = ppca_mv(dva,1,0); %missing value PCA
%     if(C(1))>0, C = -C; X = dva*C; end %first biomarker negative slope
    X = reshape(X,fliplr(size(ages)))'; X(X==0) = nan;
    X = nanmean(X,2);
    s_params = ones(numSubjects,2);
    s_params(:,1) = X; % initialize with slope = 1 and intercept = position
elseif isempty(opts.s_init) && strcmp(model,'linear')
    s_params = [zeros(size(data,1),1) ones(size(data,1),1)];
else
    s_params = opts.s_init;
    if ageoffset
        s_params(:,1) = s_params(:,1) + s_params(:,2).*age_mean;
    end
end

%--- Initial ADPS(z)
adps = bsxfun(@plus,s_params(:,1),bsxfun(@times,s_params(:,2),ages_zm));
var_k = ones(numBiomarkers,1);

%--- Initial model parameters for optimization purposes 
m_params = initializeModel(data,adps,model,bounds);
s_m = size(m_params);

%--- Initialize plot stuff
if opts.doPlots
    % subplots 2x3, 2x4 2x5 2x6 3x6 4x6...
    if numBiomarkers < 13
        ppc = 3; ppr = ceil(numBiomarkers/ppc);
%         ppc = 2; ppr = ceil(numBiomarkers/ppc);
    else
        ppr = 6; ppc = ceil(numBiomarkers/ppr);
    end
%     ppc = ppc + 2; %also plot mse convergence, subject parameters
    figure
end

% If linear model, precalculate X matrix for regression
if strcmp(model,'linear')
    Xt = zeros(numVisits*numBiomarkers,2,numSubjects);
    for i = 1:numSubjects
        % array of ages
        ages_i = ages_zm(i,:)';
        T = ones(size(ages_i,1),2);
        T(:,2) = ages_i;   
        for j = 1:numBiomarkers
            Xt((j-1)*numVisits+1:j*numVisits,:,i) = T;
        end
    end
end

%--- Optimization variables
iter = 0; err = 1e5; err_old = 1e6; sp_old = 1e5*ones(size(s_params)); mp_old = 1e5*ones(size(m_params));
% options = optimset('Display','final','TolFun',1e-6,'TolX',1e-6); 

if opts.pos_slope
    options = optimset('Display','off','Algorithm','trust-region-reflective');
    lb = [-inf 0];
else
    options = optimset('Display','off','Algorithm','levenberg-marquardt');
    lb = [];
end


% options = optimset('Display','off'); 
options2 = optimset('Display','off','LargeScale','off');
options3 = optimset('Display','off','Algorithm','active-set');


%% --- Run iterative optimization ---
while iter < opts.MaxIter
    %--- Check for convergence of parameters
%     if norm(sp_old-s_params) < opts.TolX*(eps+norm(s_params));
    if norm(mp_old-m_params) < opts.TolX*(eps+norm(m_params));
        disp(['Converged in ' num2str(iter) ' iterations'])
        break
    end
    sp_old = s_params;
    mp_old = m_params;
    fiterr = 0; err = 0;

    %% --- Fit model function to each set of biomarker data ---
    m_params = initializeModel(data,adps,model,bounds);
    for i = 1:numBiomarkers
        di = data(:,i,:);
        nv = isfinite(di(:)) & isfinite(adps(:));
        
        gamma_k_X_t = 0;
        if useCovariates
%             gamma_k = cov_params(i,:);
            gamma_k = cov_params;
            gamma_k_X = covariates*gamma_k'; % nx1
            gamma_k_X_t = repmat(gamma_k_X,[1 size(adps,2)]);
        end        
        adpsi = adps + gamma_k_X_t;
        
        if strcmp(model,'linear')
            [m_params(i,:),~,r] = regress(di(nv),[ones(sum(nv),1) adps(nv)]);
            erri = r'*r;
        else
            if isempty(bounds) || any(isnan(bounds(i,:)))
                [m_params(i,:), erri, ~, ef] = lsqcurvefit(modelfun,m_params(i,:),adpsi(nv),di(nv),[],[],options);
%                 [m_params(i,:),r,~,~,mse] = nlinfit(adps(:),di(:),modelfun,m_params(i,:));
%                 erri = nansum(r.^2);

            else    % bounds on sigmoid        
                [m_pi, erri, ~, ef] = lsqcurvefit(@(pin,xvals) sigmoidfun_bd(pin,xvals,m_params(i,[1 4])),m_params(i,2:3),adps(nv),di(nv),[],[],options);
                m_params(i,2:3) = m_pi;
            end
        end
        % Get variance of fit for weighting
        if opts.useVar
            var_k(i) = erri/(sum(nv)-s_m(2));
        end
        fiterr = 1/var_k(i)*erri + fiterr;
%         err_m = erri + err_m;
    end
    
    %--- Plot stuff
    if opts.doPlots
        % Fix X-axis bounds to remove outliers from visualization
        adps_m = nanmean(adps,2);
        qv = quantile(adps_m,[0.025 0.975]);
        % Add a buffer to the inter-quantile range
        qvr = 1*(qv(2)-qv(1));        
        qv(1) = qv(1) - qvr; qv(2) = qv(2) + qvr;
        
        for i = 1:numBiomarkers
            gamma_k_X_t = 0;
            if useCovariates
%                 gamma_k = cov_params(i,:);
                gamma_k = cov_params;
                gamma_k_X = covariates*gamma_k'; % nx1
                gamma_k_X_t = repmat(gamma_k_X,[1 size(adps,2)]);
            end        
            adpsi = adps + gamma_k_X_t;
            
            di = data(:,i,:);
            subplot(ppc,ppr,i)
            plot(adpsi(:),di(:),'c.')
            dif = feval(modelfun,m_params(i,:),adpsi(:));
            hold on, plot(adpsi(:),dif,'k.'), hold off
            title(opts.data_labels{i},'Interpreter','none','fontweight','b')
            looseAxis;
            
            % Change x-axis bounds
            a = axis;
            a(1:2) = qv;
            axis(a);
        end
        drawnow
    end
    
    %% --- Fit each subject to the model curves ---
    % Linear fit is separate since it uses simple regression
    if strcmp(model,'linear')
        % Need to add this to the main loop
        [s_params err] = fit_to_linear_model_cov(data,Xt,[],m_params,var_k); 
    else
        wts = 1./repmat(sqrt(var_k),[1 numVisits]);
%         s_params = zeros(size(s_params)); s_params(:,2) = 1;
        gam_X = zeros(numSubjects,numBiomarkers);
        if useCovariates
            gam_X = covariates*cov_params';
        end
        
        % CHANGE THIS TO PARFOR FOR PARALLEL PROCESSING
%         for i = 1:numSubjects
        parfor i = 1:numSubjects
            m_pi = m_params;
            m_pi(:,3) = m_pi(:,3) + gam_X(i,:)';
            % Get subject data
            dvi = squeeze(data(i,:,:));
            agesi = ages_zm(i,:);
            % Weight subject data
            dvi = dvi.*wts; nv = isnan(dvi); dvi(nv) = [];
            
            % Fit to model
            if sum(isfinite(agesi))>1
                [s_params(i,:), erri, ~, ef] = lsqcurvefit(@(abi,xdata) patientfun(abi,xdata,m_pi,modelfun,nv,wts),...
                            s_params(i,:),agesi,dvi',lb,[],options);
            else
                [s_params(i,:), erri, ~, ef] = lsqcurvefit(@(abi,xdata) patientfun_s(abi,xdata,m_pi,modelfun,nv,wts),...
                            s_params(i,:),agesi,dvi',[],[],options);
            end
            % Error of fit
            err = err + erri;
        end
    end    
    % Re-calculate ADPS
    adps = bsxfun(@plus,s_params(:,1),bsxfun(@times,s_params(:,2),ages_zm));
    
    %% --- Fit covariates if necessary ---
    if useCovariates
        wts = 1./sqrt(var_k);
        dc = data;
        for i = 1:numBiomarkers
            dc(:,i,:) = wts(i)*dc(:,i,:);
        end
        nv = isfinite(dc(:));
        dc = dc(nv);
        [cov_params, erri, ~, ef] = lsqcurvefit(@(covp,covs) covfun(covp,covs,adps,m_params,modelfun,wts,nv),...
                        cov_params,covariates,dc,[],[],options);
%         [cov_params, erri, ~, ef] = lsqnonlin(@(covp) covfun2(covp,covariates,adps,m_params,modelfun,wts,nv,dc),...
%             cov_params,[],[],options);
%         [cov_params, erri, ~, ef] = fminsearch(@(covp) covfun2(covp,covariates,adps,m_params,modelfun,wts,nv,dc),...
%             cov_params);
    end
    
    iter = iter+1;

    %--- Plot things again
    if opts.doPlots
%         % Plot beta vs alpha parameters of each subject
%         subplot(ppc,ppr,[(ppr*(ppc-2)+1):((ppr*(ppc-2)+1)+floor(ppr/2)-1)])
%         plot(s_params(:,1),s_params(:,2),'*')
%         xlabel('intercept (\beta)'),ylabel('slope (\alpha)')
        
        % Figure title
        set(gcf,'name',['Iter ' num2str(iter)]);
    end
    
    % Convergence values
    conv_vals = [norm(sp_old-s_params)/norm(s_params); 
                  norm(mp_old-m_params)/norm(m_params);
                  norm(fiterr-err)/norm(fiterr)];    
end
conv_vals = [iter conv_vals'];

if ageoffset
    % Correct for age offset
    s_params(:,1) = s_params(:,1) - s_params(:,2).*age_mean;
end


%% --- Function to get model ---
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

%% --- Function to fit covariates to data
function yvals = covfun(covparams,covs,adps,m_params,modelfun,wts,nv)

yvals = zeros(size(adps,1),size(m_params,1),size(adps,2));
% covparams = reshape(covparams,size(m_params,1),[]);

for i = 1:size(m_params,1)
%     gamma_k = covparams(i,:);
%     gamma_k = covparams;
%     gamma_k_X = covs*gamma_k'; % nx1
%     gamma_k_X_t = repmat(gamma_k_X,[1 size(adps,2)]);
    
    mpi = m_params(i,:); mpi(3) = mpi(3) + covparams;
    yvals(:,i,:) = wts(i).*feval(modelfun,mpi,adps);  
%     yvals(:,i,:) = wts(i).*feval(modelfun,m_params(i,:),adps+gamma_k_X_t);    
end
% nv = isfinite(yvals(:));
yvals = yvals(nv);

function yvals = covfun2(covparams,covs,adps,m_params,modelfun,wts,nv,data)

yvals = zeros(size(adps,1),size(m_params,1),size(adps,2));
% covparams = reshape(covparams,size(m_params,1),[]);

for i = 1:size(m_params,1)
%     gamma_k = covparams(i,:);
    gamma_k = covparams;
    gamma_k_X = covs*gamma_k'; % nx1
    gamma_k_X_t = repmat(gamma_k_X,[1 size(adps,2)]);
    
    yvals(:,i,:) = wts(i).*feval(modelfun,mpi,adps+gamma_k_X_t);    
end
% nv = isfinite(yvals(:));
yvals = yvals(nv);
yvals = yvals-data;
yvals = yvals'*yvals;


%% --- Function to evaluate fit of subject to model ---
function yvals = patientfun(a_coefs,ages,m_params,modelfun,fv,wts)
% Note jacobian would be B(y-D)*(1-(y-D)/A) and t*B(y-D)*(1-(y-D)/A)

adps = a_coefs(1) + ages*a_coefs(2);
yvals = feval(modelfun,m_params,adps);
yvals = yvals.*wts;
yvals(fv) = [];
yvals = yvals';


%% --- Function to evaluate fit of subject to model - single time point ---
function yvals = patientfun_s(a_coefs,ages,m_params,modelfun,fv,wts)
% Note jacobian would be B(y-D)*(1-(y-D)/A) and t*B(y-D)*(1-(y-D)/A)

adps = a_coefs(1) + ages;
yvals = feval(modelfun,m_params,adps);
yvals = yvals.*wts;
yvals(fv) = [];
yvals = yvals';


%% --- Function to initialize model fit ---
function m_params = initializeModel(data,adps,model,bounds)

[~, numBiomarkers, ~] = size(data);
switch model
    case 'linear'
        %--- Linear model, dont need to initialize
        m_params = zeros(numBiomarkers,2);
    case 'sigmoid'
        %--- Sigmoid fit, try to be smart about center, max and min based
        %    on current ADPS estimate
        m_params = ones(numBiomarkers,4);
        for i = 1:numBiomarkers
            % biomarker data
            yvi = squeeze(data(:,i,:));
            zvi = adps(isfinite(yvi));
            yvi = yvi(isfinite(yvi));
            
            % max and min are obvious...
            m_params(i,4) = nanmin(yvi(:));
            m_params(i,1) = nanmax(yvi(:))-m_params(i,4);
            
            % sort by ADPS and make slope negative if necessary
            [~, indzs] = sort(zvi);
            if nanmean(yvi(indzs(1:20))) > nanmean(yvi(indzs(end-20:end)))
                m_params(i,2) = -m_params(i,2);
            end
            
            % cut data into 'positive' and 'negative' values based on
            %   halfway point of max and min values
            z_p = nanmedian(zvi(yvi>(nanmax(yvi(:))+nanmin(yvi(:)))/2));
            z_n = nanmedian(zvi(yvi<(nanmax(yvi(:))+nanmin(yvi(:)))/2));
            m_params(i,3) = -mean([z_p,z_n]);
            
            % if there are bounds, just use those
            if ~isempty(bounds)
                if ~any(isnan(bounds(i,:)))
                    m_params(i,1) = max(bounds(i,:))-min(bounds(i,:));
                    m_params(i,4) = min(bounds(i,:));
                end
            end
        end
    case 'richards'
        % Initialize same as sigmoid
        m_params = ones(numBiomarkers,5);
        m_params(:,1:4) = initializeModel(data,adps,'sigmoid',[]);
    case 'quadratic'
        m_params = ones(numBiomarkers,3);
    otherwise
        error(['Invalid model: ' model]);
end


%% --- Get options and set empty to defaults ---
function newopts = getopts(opts)

% Defaults
options = struct('MaxIter', 100, 'TolFun', 1e-7, 'TolX', 1e-7,...
    'useVar', false, 'ar1model', false, 'doPlots', false,'s_init',[],...
    'bounds',[],'ageOffset',false,'model','linear','covariates',[],...
    'data_labels',[],'pos_slope',false);

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


%% --- Fit function for bounded sigmoid model ---
function yvals = sigmoidfun_bd(pin,xvals,maxmin)
% put in new file?
% bounded sigmoid function, 
%   y = a/(1 + exp(-b(x+c))) + d with a and d fixed

if any(size(pin) == 1)
    a = maxmin(1); d = maxmin(2);
    yvals = a./(1+exp(-pin(1)*(xvals+pin(2)))) + d;
else
    a = sum(maxmin,2); d = min(maxmin,[],2);
    yvals = bsxfun(@plus,bsxfun(@rdivide,a,(1+exp(bsxfun(@times,-pin(:,1),bsxfun(@plus,xvals,pin(:,2)))))),d);
end

