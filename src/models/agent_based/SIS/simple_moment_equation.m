function dydt = simple_moment_equation(t,y,system)
% The constraint equation in forward form

% compute meanfield infection rate
phi_SI = system.rates(1)*(1-y).*(system.graph*y);

% initialize
dydt = -system.rates(2)*y+phi_SI;



end

