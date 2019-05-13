function dydt = forward_equation(t,m,alpha_t,alpha)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(9,1);
dydt(1) = alpha(1)*(1-m(1))-alpha(2)*m(1);
dydt(2) = alpha(3)*m(1)-alpha(4)*m(2);
dydt(3) = alpha(5)*m(2)-alpha(6)*m(3);
dydt(4) = -alpha(1)*m(1)+alpha(1)-2*alpha(1)*m(4)-2*alpha(2)*m(4)+alpha(2)*m(1);
dydt(5) = -alpha(1)*m(5)-alpha(2)*m(5)+alpha(3)*m(4)-alpha(4)*m(5);
dydt(6) = -alpha(1)*m(6)-alpha(2)*m(6)+alpha(5)*m(5)-alpha(6)*m(6);
dydt(7) = 2*alpha(3)*m(5)+alpha(3)*m(1)-2*alpha(4)*m(7)+alpha(4)*m(2);
dydt(8) = alpha(3)*m(6)-alpha(4)*m(8)+alpha(5)*m(7)-alpha(6)*m(8);
dydt(9) = 2*alpha(5)*m(8)+alpha(5)*m(2)-2*alpha(6)*m(9)+alpha(6)*m(3);

end

