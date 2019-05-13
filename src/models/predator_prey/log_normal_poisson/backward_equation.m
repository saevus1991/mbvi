function dydt = backward_equation(t,eta,alpha_t,alpha,m,c)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(5,1);
dydt(1) = c(1,1)-alpha(1)+alpha(1)*log(alpha(1)/c(1,2))-alpha(3)*eta(2)*m(2)+(c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2)))*m(2)+(c(3,1)-alpha(3)+alpha(3)*log(alpha(3)/c(3,2)))*m(2)-eta(1)*(alpha(1)-alpha(2)*m(2))-eta(3)*(alpha(1)+alpha(2)*m(2)-2*alpha(2)*(m(2)-2*m(1)*m(2)-m(4)+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/m(1)^2+((-1+2*m(1))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2))-(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^3*m(2))))-eta(4)*(-alpha(3)*m(1)*m(2)+alpha(2)*m(2)^2-alpha(3)*(m(1)*m(2)+m(4))+alpha(3)*(m(2)+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/m(1)^2+((-1+2*m(1))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2))-(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^3*m(2)))-alpha(2)*(m(2)+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2))-((m(1)*m(2)+m(4))^2*(-m(2)+m(2)^2+m(5)))/(m(1)^2*m(2)^2)))-eta(5)*(alpha(3)*m(2)+2*alpha(3)*(m(2)-m(2)^2+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2))-((m(1)*m(2)+m(4))^2*(-m(2)+m(2)^2+m(5)))/(m(1)^2*m(2)^2)));
dydt(2) = c(4,1)-alpha(4)+alpha(4)*log(alpha(4)/c(4,2))+alpha(2)*eta(1)*m(1)+(c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2)))*m(1)+(c(3,1)-alpha(3)+alpha(3)*log(alpha(3)/c(3,2)))*m(1)-eta(2)*(-alpha(4)+alpha(3)*m(1))-eta(3)*(alpha(2)*m(1)-2*alpha(2)*(m(1)-m(1)^2+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)*m(2))-((-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)^2)))-eta(4)*(-alpha(3)*m(1)^2+alpha(2)*m(1)*m(2)+alpha(2)*(m(1)*m(2)+m(4))+alpha(3)*(m(1)+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)*m(2))-((-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)^2))-alpha(2)*(m(1)+((-1+2*m(2))*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2)+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/m(2)^2-(2*(m(1)*m(2)+m(4))^2*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2)^3)))-eta(5)*(alpha(4)+alpha(3)*m(1)+2*alpha(3)*(m(1)-2*m(1)*m(2)-m(4)+((-1+2*m(2))*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2)+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/m(2)^2-(2*(m(1)*m(2)+m(4))^2*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2)^3)));
dydt(3) = -((alpha(3)*eta(4)*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)))-eta(3)*(2*alpha(1)-(2*alpha(2)*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)));
dydt(4) = c(2,1)+c(3,1)-alpha(2)-alpha(3)+alpha(2)*eta(1)-alpha(3)*eta(2)+alpha(2)*log(alpha(2)/c(2,2))+alpha(3)*log(alpha(3)/c(3,2))-eta(3)*(alpha(2)-2*alpha(2)*(1-m(1)+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)^2*m(2))))-eta(4)*(alpha(1)-alpha(4)-alpha(3)*m(1)+alpha(2)*m(2)+alpha(3)*(1+(2*(-m(1)+m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)^2*m(2)))-alpha(2)*(1+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2)^2)))-eta(5)*(alpha(3)+2*alpha(3)*(1-m(2)+(2*(m(1)*m(2)+m(4))*(-m(2)+m(2)^2+m(5)))/(m(1)*m(2)^2)));
dydt(5) = (alpha(2)*eta(4)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2)-eta(5)*(-2*alpha(4)+(2*alpha(3)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2));

end