function dydt = forward_gradient_evolution(t,m,c)
% An ode for the prior dynamics for fixed parameters along with the
% jacobian matrix with respect to the parameters

% evaluate the derivative
dydt = zeros(9,1);
dydt(1) = c(1)*(1-m(1))-c(2)*m(1);
dydt(2) = c(3)*m(1)-c(4)*m(2);
dydt(3) = c(5)*m(2)-c(6)*m(3);

% evaluate the jacobian elements

% dm1/dc1
dydt(4) = (1-m(1))-(c(1)+c(2))*m(4);
% dm1/dc2
dydt(5) = -m(1)-(c(1)+c(2))*m(5);
% dm1/dc3
dydt(6) = 0;
% dm1/dc4 
dydt(7) = 0;
% dm1/dc5
dydt(8) = 0;
% dm1/dc6
dydt(9) = 0;

% dm2/dc1 
dydt(10) = c(3)*m(4)-c(4)*m(10);
% dm2/dc2
dydt(11) = c(3)*m(5)-c(4)*m(11);
% dm2/dc3
dydt(12) = m(1)+c(3)*0-c(4)*m(12);
% dm2/dc4
dydt(13) = -m(2)+c(3)*0-c(4)*m(13);
% dm2/dc5
dydt(14) = 0;
% dm2/dc6
dydt(15) = 0;

% dm3/dc1
dydt(16) = c(5)*m(10)-c(6)*m(16);
% dm3/dc2
dydt(17) = c(5)*m(11)-c(6)*m(17);
% dm3/dc3
dydt(18) = c(5)*m(12)-c(6)*m(18);
% dm3/dc4
dydt(19) = c(5)*m(13)-c(6)*m(19);
% dm3/dc5
dydt(20) = m(2)+c(5)*0-c(6)*m(20);
% dm3/dc6
dydt(21) = -m(3)+c(5)*0-c(6)*m(21);

end

