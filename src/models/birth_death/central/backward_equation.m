function dydt = backward_equation(t,eta,alpha_t,alpha,c)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(2,1);
dydt(1) = c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2))+alpha(2)*eta(1)...
    -alpha(2)*eta(2);
dydt(2) = 2*alpha(2)*eta(2);

end

