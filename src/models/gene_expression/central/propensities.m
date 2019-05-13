function [prop] = propensities(m)
% The constraint equation in forward form

% set up output
prop = zeros(size(m,1),6);

% gene on
prop(:,1) = (1-m(:,1));

% gene off
prop(:,2) = m(:,1);

% translation
prop(:,3) = m(:,1);

% mrna degradation
prop(:,4) = m(:,2);

% transcription
prop(:,5) = m(:,2);

% protein degradation
prop(:,6) = m(:,3);


end

