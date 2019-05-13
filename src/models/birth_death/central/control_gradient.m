function [grad] = control_gradient(eta,m)
% The constraint equation in forward form

% evalute the contribution to the gradient
grad = zeros(size(m,1),2);
grad(:,1) = eta(:,1)+eta(:,2);
grad(:,2) = -m(:,1).*eta(:,1)-2*m(:,2).*eta(:,2)+eta(:,2).*m(:,1);

end

