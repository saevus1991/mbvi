function [grad] = control_gradient(eta,m)
% The constraint equation in forward form

% evalute the contribution to the gradient
grad = zeros(size(m,1),6);
grad(:,1) = -eta(:,1).*(1-m(:,1))-eta(:,4).*(1-m(:,1)-2.*m(:,4))+eta(:,5).*m(:,5)+eta(:,6).*m(:,6);
grad(:,2) = eta(:,1).*m(:,1)-eta(:,4).*(m(:,1)-2.*m(:,4))+eta(:,5).*m(:,5)+eta(:,6).*m(:,6);
grad(:,3) = -eta(:,2).*m(:,1)-eta(:,5).*m(:,4)-eta(:,7).*(m(:,1)+2.*m(:,5))-eta(:,8).*m(:,6);
grad(:,4) = eta(:,2).*m(:,2)+eta(:,5).*m(:,5)-eta(:,7).*(m(:,2)-2.*m(:,7))+eta(:,8).*m(:,8);
grad(:,5) = -eta(:,3).*m(:,2)-eta(:,6).*m(:,5)-eta(:,8).*m(:,7)-eta(:,9).*(m(:,2)+2.*m(:,8));
grad(:,6) = eta(:,3).*m(:,3)+eta(:,6).*m(:,6)+eta(:,8).*m(:,8)-eta(:,9).*(m(:,3)-2.*m(:,9));

grad = -grad;

end