function dydt = forward_equation(t,m,alpha_t,alpha)
% The constraint equation in forward form

% perform interpolation
alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));

% evaluate the derivative
dydt = zeros(5,1);
dydt(1) = alpha(1)*m(1)-alpha(2)*(m(4)+m(1)*m(2));
dydt(2) = alpha(3)*(m(4)+m(1)*m(2))-alpha(4)*m(2);
dydt(3) = 2*alpha(1)*m(3)+alpha(1)*m(1)+alpha(2)*(m(4)+m(1)*m(2))-2*alpha(2)*(m112(m)-(m(4)+m(1)*m(2))*m(1));
dydt(4) = alpha(1)*m(4)-alpha(2)*m122(m)+alpha(2)*(m(4)+m(1)*m(2))*m(2)+alpha(3)*m112(m)-alpha(3)*(m(4)+m(1)*m(2))*m(1)-alpha(4)*m(4);
dydt(5) = 2*alpha(3)*(m122(m)-(m(4)+m(1)*m(2))*m(2))+alpha(3)*(m(4)+m(1)*m(2))-2*alpha(4)*m(5)+alpha(4)*m(2);

end

function out = m112(m)

out = ((m(3)+m(1)*m(1)-m(1))/m(2))*((m(4)+m(1)*m(2))/m(1))^2+m(4)+m(1)*m(2);

end

function out = m122(m)

out = ((m(5)+m(2)*m(2)-m(2))/m(1))*((m(4)+m(1)*m(2))/m(2))^2+m(4)+m(1)*m(2);

end