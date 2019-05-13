%% test ssa for SIS disease model on a graph

addpath('../src/ssa')

% specify graph by adjecency matrix

A = [0,1,0,0,0,0;
     1,0,1,1,0,0;
     0,1,0,0,0,0;
     0,1,0,0,1,1;
     0,0,0,1,0,0;
     0,0,0,1,0,0];
 
 % initial state
 x = [1;0;0;0;0;0];
 rates = [0.5;0.1];
 
 % set up a system structure for stochastic simulation
 system.initial = x;
 system.rates = rates;      % infection rate, recovery rate
 system.graph = A;
 system.dynamcis = 'SIS';
 system.t_min = 0;
 system.t_max = 100;
 
 % simulate example trajectory
 
 [times,states] = ssa_sis(system);
 
