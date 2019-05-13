function [prop] = propensities(m)
% The constraint equation in forward form

% set up output
prop = zeros(size(m,1),4);

% gene on
prop(:,1) = m(:,1);

% gene off
prop(:,2) = m(:,4)+m(:,1).*m(:,2);

% translation
prop(:,3) = m(:,4)+m(:,1).*m(:,2);

% mrna degradation
prop(:,4) = m(:,2);

end

