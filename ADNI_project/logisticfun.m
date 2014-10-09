function [yvals jacob] = logisticfun(pin,xvals)
% function yvals = logisticfun(pin,xvals)
% y = a/(1 + exp(-b(x+c))) + d
if any(size(pin) == 1)
    yvals = pin(1)./(1+exp(-pin(2)*(xvals+pin(3)))) + pin(4);
else
    yvals = bsxfun(@plus,bsxfun(@rdivide,pin(:,1),(1+exp(bsxfun(@times,-pin(:,2),bsxfun(@plus,xvals,pin(:,3)))))),pin(:,4));
end

if nargout > 1
    jacob = zeros(length(pin),length(xvals));
    jacob(1,:) = 1./(1+exp(-pin(2)*(xvals+pin(3))));
    expv = exp(-pin(2)*(xvals+pin(3)))./(1+exp(-pin(2)*(xvals+pin(3)))).^2;
    if any(isnan(expv))
        expv(isnan(expv)) = 0;
    end
    jacob(2,:) = (xvals+pin(3))*pin(1).*expv;
    jacob(3,:) = pin(2)*pin(1)*expv;
    jacob(4,:) = 1;
    jacob = jacob';
end