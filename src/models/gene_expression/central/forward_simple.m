function dydt = forward_simple(t,m,c)
% The prio dynamics just for fixed parameters

% evaluate the derivative
dydt = zeros(9,1);
dydt(1) = c(1)*(1-m(1))-c(2)*m(1);
dydt(2) = c(3)*m(1)-c(4)*m(2);
dydt(3) = c(5)*m(2)-c(6)*m(3);
dydt(4) = -c(1)*m(1)+c(1)-2*c(1)*m(4)-2*c(2)*m(4)+c(2)*m(1);
dydt(5) = -c(1)*m(5)-c(2)*m(5)+c(3)*m(4)-c(4)*m(5);
dydt(6) = -c(1)*m(6)-c(2)*m(6)+c(5)*m(5)-c(6)*m(6);
dydt(7) = 2*c(3)*m(5)+c(3)*m(1)-2*c(4)*m(7)+c(4)*m(2);
dydt(8) = c(3)*m(6)-c(4)*m(8)+c(5)*m(7)-c(6)*m(8);
dydt(9) = 2*c(5)*m(8)+c(5)*m(2)-2*c(6)*m(9)+c(6)*m(3);

end

