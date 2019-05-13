% load required libraries
addpath('../src/ssa')
addpath('../src/models/predator_prey')
%addpath('models/predator_prey/log_normal')

%% Description
% Compare the reaction rate equations of the predator prey model with
% simulation results and state space truncation

%% Initialize parameters

% construct priors
rates = [5e-4;1e-4;1e-4;5e-4];

% observation parameters
sigma = 10;            % observation noise (gaussian mode

%% construct the model

Pre = [1,0; ...
      1,1; ...
      1,1; ...
      0,1];
Post = [2,0; ...
       0,1; ...
       1,2; ...
       0,0];

system = struct;
system.name = 'Simple Eucaryotic Gene Expression Model';
system.Pre = Pre;
system.Post = Post;
system.S = Post-Pre;
system.rates = rates;
system.initial = [20;5];
system.max_state = [100;100];
system.t_min = 0;
system.t_max = 3000;
system.num_species = size(system.Pre,2);
system.num_reactions = size(system.Pre,1);
system.seed = randi(intmax);

% Gaussian observation noise
noise = @(x) normrnd(x,sigma);

%% simulation

% draw a sample path
[reaction_time,state] = gillespie(system);

% discretize the jump process trajectories for plotting
delta = 1;
N = 10000;
t_discrete = [system.t_min:delta:system.t_max-delta];
moments = zeros(4,length(t_discrete));
for i = 1:N
   system.seed = randi(intmax);
   [reaction_time,state] = gillespie(system);
   [discrete_state] = discretize_data(reaction_time,state,t_discrete); 
   moments(1:2,:) = moments(1:2,:)+discrete_state;
   moments(3:4,:) = moments(1:2,:)+discrete_state.^2;
end
moments = moments/N;

% alpha_t = [system.t_min;system.t_max];
% alpha = [rates';rates'];
% moment_ode = @(t,y) forward_equation(t,y,alpha_t,alpha);
%  [sol_t,sol_state] = ode45(moment_ode,t_discrete,[system.initial;0;0;0]);
 
rate_equation = @(t,m) [rates(1)*m(1)-rates(2)*m(1)*m(2);rates(3)*m(1)*m(2)-rates(4)*m(2)];
 [rre_t,rre_state] = ode45(rate_equation,t_discrete,system.initial);

 %% integration with state space truncation
 
 % set up integrator object
 generator = Generator_PP(system);
 
 % get initial distribution
 initial_dist = generator.initial_dist(system.initial);
 
 % solve master equation
 [t,dist] = generator.integrate(initial_dist,system.t_min,system.t_max);
 
 % compute marginal
 population_mean = generator.population_mean(dist);

 %% plotting
 
 % simulated mean vs rre mean
 figure
 hold on
 plot(t_discrete,moments(1,:),'-b','linewidth',2)
 plot(rre_t,rre_state(:,1),'--b','linewidth',2)
 plot(t_discrete,moments(2,:),'-r','linewidth',2)
 plot(rre_t,rre_state(:,2),'--r','linewidth',2)
 ylim([0,50])
 
 % simulated mean vs truncated mean
 figure
 hold on
 plot(t_discrete,moments(1,:),'-b','linewidth',2)
 plot(t,population_mean(1,:),'--b','linewidth',2)
 plot(t_discrete,moments(2,:),'-r','linewidth',2)
 plot(t,population_mean(2,:),'--r','linewidth',2)
 ylim([0,50])
