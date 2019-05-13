%% test ssa for SIS disease model on a graph

addpath('../src/models/autocatalytic_birth_death')
addpath('../src/ssa')
addpath('../src/mean_field_vi')

%% set up predator prey model
 
 % set up a system structure for stochastic simulation
 pp_system.initial = [10;5];
 pp_system.max_state = [100;100];
 pp_system.rates = [5e-4;1e-4;1e-4;5e-4];      % infection rate, recovery rate
 pp_system.t_min = 0;
 pp_system.t_max = 10000;
 
 % construct time grid
 delta_t = 0.5;
 time = pp_system.t_min:delta_t:pp_system.t_max;
 time(end) = time(end)-0.001*delta_t;
 
 %% initialize forward parts of the components
 
 % set up autocatalytic process as template for prey
 prey_system.initial = pp_system.initial(1);
 prey_system.max_state = pp_system.max_state(1);
 prey_system.rates = pp_system.rates(1:2);      % infection rate, recovery rate
 prey_system.t_min = pp_system.t_min;
 prey_system.t_max = pp_system.t_max;
 
 % set up autocatalytic process as template for predator
 pred_system.initial = pp_system.initial(2);
 pred_system.max_state = pp_system.max_state(2);
 pred_system.rates = pp_system.rates(3:4);      % infection rate, recovery rate
 pred_system.t_min = pp_system.t_min;
 pred_system.t_max = pp_system.t_max;
 
 % construct static generators
 prey_generator = Generator_ABD(prey_system);
 pred_generator = Generator_ABD(pred_system);
 
% construct time dependent forward generators
prey_forward = Generator(time,prey_generator.generator);
pred_forward = Generator(time,pred_generator.generator);

% test initialization
initial_dist = prey_generator.initial_dist(prey_system.initial);
[prey_dist,~] = prey_forward.integrate(initial_dist);
initial_dist = pred_generator.initial_dist(pred_system.initial);
[pred_dist,~] = pred_forward.integrate(initial_dist);
prey_mean = prey_generator.population_mean(prey_dist);
pred_mean = pred_generator.population_mean(pred_dist);
 
% plot means
 % plot marginals
 figure
 hold on
 plot(time,prey_mean,'-r')
 plot(time,pred_mean,'-b')
 xlabel('Time','fontsize',14)

 % set up a system structure for stochastic simulation

 
%  % set up integrator object
%  generator = Generator_ABD(system);
%  
%  % get initial distribution
%  initial_dist = generator.initial_dist(system.initial);
%  
%  % solve master equation
%  [t,dist] = generator.integrate(initial_dist,system.t_min,system.t_max);
%  
%  % compute marginal
%  population_mean = generator.population_mean(dist);
%  
%  % plot marginals
%  figure
%  hold on
%  plot(t,population_mean,'-r')
%  xlabel('Time','fontsize',14)
%  
%  % simulate example trajectory
%  
%  %[times,states] = ssa_sis(system);
%  
