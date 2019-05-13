function dydt = forward_equation(t,m,alpha_t,alpha)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(2,1);
dydt(1) = alpha(1)-alpha(2)*m(1);
dydt(2) = alpha(1)-2*alpha(2)*m(2)+alpha(2)*m(1);

end

