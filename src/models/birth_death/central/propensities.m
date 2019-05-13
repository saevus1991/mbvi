function [prop] = propensities(m)
% The constraint equation in forward form

% set up output
prop = zeros(size(m,1),2);

% birth reaction
prop(:,1) = 1;

% death reaction
prop(:,2) = m(:,1);

end

