%% test ssa for SIS disease model on a graph

addpath('../src/models/predator_prey')
addpath('../src/ssa')

% specify graph by adjecency matrix
 
 % initial state
 x = [10;5];
 rates = [5e-4;1e-4;1e-4;5e-4];
 max_state = [100;100];
 
 % set up a system structure for stochastic simulation
 system.initial = x;
 system.max_state = max_state;
 system.rates = rates;      % infection rate, recovery rate
 system.t_min = 0;
 system.t_max = 10000;
 
 % set up integrator object
 generator = Generator_PP(system);
 
 % get initial distribution
 initial_dist = generator.initial_dist(system.initial);
 
 % solve master equation
 [t,dist] = generator.integrate(initial_dist,system.t_min,system.t_max);
 
 % compute marginal
 population_mean = generator.population_mean(dist);
 
 % plot marginals
 figure
 hold on
 plot(t,population_mean(1,:),'-r')
 plot(t,population_mean(2,:),'-b')
 xlabel('Time','fontsize',14)
 
 % simulate example trajectory
 
 %[times,states] = ssa_sis(system);
 
