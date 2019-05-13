function [grad] = control_gradient(eta,m)
% The constraint equation in forward form

% evalute the contribution to the gradient
grad = zeros(size(m,1),4);
grad(:,1) = -eta(:,1).*m(:,1)-eta(:,3).*(m(:,1)+2.*m(:,3))-eta(:,4).*m(:,4);
grad(:,2) = -eta(:,1).*(-m(:,1).*m(:,2)-m(:,4))-eta(:,3).*(m(:,1).*m(:,2)+m(:,4)-2.*(m(:,2).*(-m(:,1).^2+m(:,3))+m(:,1).*(m(:,1).*m(:,2)+m(:,4))))-eta(:,4).*(-m(:,2).*(m(:,1).*m(:,2)+m(:,4))-m(:,1).*(-m(:,2).^2+m(:,5)));
grad(:,3) = -eta(:,2).*(m(:,1).*m(:,2)+m(:,4))-eta(:,4).*(m(:,2).*(-m(:,1).^2+m(:,3))+m(:,1).*(m(:,1).*m(:,2)+m(:,4)))-eta(:,5).*(m(:,1).*m(:,2)+m(:,4)+2.*(m(:,2).*(m(:,1).*m(:,2)+m(:,4))+m(:,1).*(-m(:,2).^2+m(:,5))));
grad(:,4) = eta(:,2).*m(:,2)+eta(:,4).*m(:,4)-eta(:,5).*(m(:,2)-2.*m(:,5));

grad = -grad;

end