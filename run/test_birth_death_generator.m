%% test ssa for SIS disease model on a graph

addpath('../src/models/birth_death')
addpath('../src/ssa')

% specify graph by adjecency matrix
 
 % initial state
 x = 0;
 rates = [0.5;0.05];
 max_state = 100;
 
 % set up a system structure for stochastic simulation
 system.initial = x;
 system.max_state = max_state;
 system.rates = rates;      % infection rate, recovery rate
 system.t_min = 0;
 system.t_max = 100;
 
 % set up integrator object
 generator = Generator_BD(system);
 
 % get initial distribution
 initial_dist = generator.initial_dist(system.initial);
 
 % solve master equation
 [t,dist] = generator.integrate(initial_dist,system.t_min,system.t_max);
 
 % compute marginal
 population_mean = generator.population_mean(dist);
 
 % plot marginals
 figure
 hold on
 plot(t,population_mean,'-r')
 xlabel('Time','fontsize',14)
 
 % simulate example trajectory
 
 %[times,states] = ssa_sis(system);
 
