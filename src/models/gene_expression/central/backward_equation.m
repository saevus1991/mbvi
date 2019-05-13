function dydt = backward_equation(t,eta,alpha_t,alpha,m,c)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(9,1);
dydt(1) = -c(1,1)+c(2,1)+c(3,1)+alpha(1)-alpha(2)-alpha(3)-(-alpha(1)-alpha(2))*eta(1)-alpha(3)*eta(2)-(-alpha(1)+alpha(2))*eta(4)-alpha(3)*eta(7)-alpha(1)*log(alpha(1)/c(1,2))+alpha(2)*log(alpha(2)/c(2,2))+alpha(3)*log(alpha(3)/c(3,2));
dydt(2) = c(4,1)+c(5,1)-alpha(4)-alpha(5)+alpha(4)*eta(2)-alpha(5)*eta(3)-alpha(4)*eta(7)-alpha(5)*eta(9)+alpha(4)*log(alpha(4)/c(4,2))+alpha(5)*log(alpha(5)/c(5,2));
dydt(3) = c(6,1)-alpha(6)+alpha(6)*eta(3)-alpha(6)*eta(9)+alpha(6)*log(alpha(6)/c(6,2));
dydt(4) = -(-2*alpha(1)-2*alpha(2))*eta(4)-alpha(3)*eta(5);
dydt(5) = -(-alpha(1)-alpha(2)-alpha(4))*eta(5)-alpha(5)*eta(6)-2*alpha(3)*eta(7);
dydt(6) = -(-alpha(1)-alpha(2)-alpha(6))*eta(6)-alpha(3)*eta(8);
dydt(7) = 2*alpha(4)*eta(7)-alpha(5)*eta(8);
dydt(8) = -(-alpha(4)-alpha(6))*eta(8)-2*alpha(5)*eta(9);
dydt(9) = 2*alpha(6)*eta(9);

end