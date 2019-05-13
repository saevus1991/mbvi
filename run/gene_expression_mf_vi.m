% load required libraries
addpath('../src/ssa')
addpath('../src/models/gene_expression')
%addpath('models/predator_prey/log_normal')

%% Description
% Apply the mean-field variational inference to the gene expression model

%% Initialize parameters

% fix rng seed for reproducibility
rng(3856820-2)

% set rates
rates = [0.001;0.001;0.15;0.001;0.04;0.008];

% observation parameters
sigma = 15;            % observation noise (gaussian mode

%% construct the model

Pre = [1,0,0,0; ...
      0,1,0,0; ...
      0,1,0,0; ...
      0,0,1,0;...
      0,0,1,0;...
      0,0,0,1];
Post = [0,1,0,0; ...
       1,0,0,0; ...
       0,1,1,0; ...
       0,0,0,0; ...
       0,0,1,1; ...
       0,0,0,0];

system = struct;
system.name = 'Simple Eucaryotic Gene Expression Model';
system.Pre = Pre;
system.Post = Post;
system.S = Post-Pre;
system.rates = rates;
system.initial = [0;1;0;0];
system.t_min = 0;
system.t_max = 5200;
system.num_species = size(system.Pre,2);
system.num_reactions = size(system.Pre,1);
system.seed = randi(intmax);

% Gaussian observation noise
noise = @(x) normrnd_2(x,sigma);

% constants
delta = 1;
delta_t = 400;


%% simulation

% % draw a sample path
% [reaction_time,state] = gillespie(system);
% 
% % discretize the jump process trajectories for plotting
% tic
% delta = 1;
% N = 10000;
% t_discrete = [system.t_min:delta:system.t_max-delta];
% moments = zeros(3,length(t_discrete));
% for i = 1:N
%    system.seed = randi(intmax);
%    [reaction_time,state] = gillespie(system);
%    [discrete_state] = discretize_data(reaction_time,state,t_discrete); 
%    moments = moments+discrete_state(2:end,:);
% end
% moments = moments/N;
% end_time = toc;
% fprintf("Simulation required %f seconds \n",end_time);
%  
% rate_equation = @(t,m) [rates(1)*(1-m(1))-rates(2)*m(1);rates(3)*m(1)-rates(4)*m(2);rates(5)*m(2)-rates(6)*m(3)];
%  [rre_t,rre_state] = ode45(rate_equation,t_discrete,system.initial(2:end));



 %% plotting
 
%  % simulated mean vs rre mean
%  figure
%  hold on
%  plot(t_discrete,moments(2,:),'-b','linewidth',2)
%  plot(rre_t,rre_state(:,2),'--b','linewidth',2)
%  plot(t_discrete,moments(3,:),'-r','linewidth',2)
%  plot(rre_t,rre_state(:,3),'--r','linewidth',2)
%  %ylim([0,50])
 
%  % simulated mean vs truncated mean
%  figure
%  hold on
%  plot(t_discrete,moments(1,:),'-b','linewidth',2)
%  plot(t,population_mean(1,:),'--b','linewidth',2)
%  plot(t_discrete,moments(2,:),'-r','linewidth',2)
%  plot(t,population_mean(2,:),'--r','linewidth',2)
%  ylim([0,50])
