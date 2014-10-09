function m_params = convertLogisticParams(m_params_in)
% We want consistancy with the logistic function formulation
% Noting that A + B/(1+exp(-C(x+D))) = A + B - B/(1+exp(C(x+D)))
% We want the first formulation, with positive C value when _/`` and
%   negative C value when ``\_

m_params = m_params_in;

% find inconsistant notations (B < 0)
inds = m_params_in(:,1) < 0;
m_params(inds,4) = m_params_in(inds,1) + m_params_in(inds,4);
m_params(inds,1) = -m_params_in(inds,1);
m_params(inds,2) = -m_params_in(inds,2);