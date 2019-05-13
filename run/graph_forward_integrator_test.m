%% test ssa for SIS disease model on a graph

addpath('../src/models/agent_based/sis')
addpath('../src/ssa')

% set up predator prey model

% specify graph by adjecency matrix

A = [0,1,0,0,0,0;
     1,0,1,1,0,0;
     0,1,0,0,0,0;
     0,1,0,0,1,1;
     0,0,0,1,0,0;
     0,0,0,1,0,0];
 
 % initial state
 x = [1;0;0;0;0;0];
 rates = [0.05;0.1];
 
 % set up a system structure for stochastic simulation
 system.initial = x;
 system.rates = rates;      % infection rate, recovery rate
 system.graph = A;
 system.dynamcis = 'SIS';
 system.t_min = 0;
 system.t_max = 100;
 
 % set up integrator object 
 integrator = Generator_SIS(system);
 
 % get initial distribution
 initial_dist = integrator.initial_dist(system.initial);

 % solve master equation
 [t,dist] = integrator.integrate(initial_dist,system.t_min,system.t_max);
 
 % compute marginal
 marginal_dist = integrator.marginalize(dist);
 
 % compute meanfield solutions of the marginals
 ode_fun = @(t,m) simple_moment_equation(t,m,system);
 [t_mf,marginal_mf] = ode45(ode_fun,[system.t_min,system.t_max],system.initial);
 marginal_mf = marginal_mf';
 
 % compute a sample average using simulations
 N = 1000;
 marginal_sim = zeros(length(system.initial),N);
 t_sample = linspace(system.t_min,system.t_max,1000);
 for i = 1:N
    [times,states] = ssa_sis(system);
    discrete_data = discretize_data([system.t_min,times],[system.initial,states],t_sample);
    marginal_sim = marginal_sim+discrete_data;
 end
 marginal_sim = marginal_sim/N;
 
 % plot marginals
 figure
 for i = 1:6
 subplot(3,2,i)
 hold on
 plot(t,marginal_dist(i,:),'-r')
 plot(t_mf,marginal_mf(i,:),'-b')
  plot(t_sample,marginal_sim(i,:),'-g')
 ylabel(num2str(i),'fontsize',14)
 xlabel('Time','fontsize',14)
 end
 
 
 % simulate example trajectory
 
 %[times,states] = ssa_sis(system);
 
