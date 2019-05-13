% load required libraries
addpath('/Users/christian/Documents/Latex/Publications/ICML2019/post_submission/code/src/ssa')

%% Description
% Simulate N trajectories from the gene exrpession model with log-normal
% noise

%% Initialize parameters

% fix random number seed for reproducibility
rng(6868183)

% set rates
rates = [0.001;0.001;0.15;0.001;0.04;0.008];

% number of trajectories
N = 100;

%% construct the model

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
sigma = 0.15;  
noise = @(x) lognrnd(log(x),sigma);

% constants
delta = 1;
delta_t = 400;

%% simulation

% observed species
selection = [4];

% initilize dataset
data_set.system = system;
data_set.noise_model = noise;
data_set.sigma = sigma;
data_set.reaction_times = cell(N,1);
data_set.states = cell(N,1);
data_set.obs_times = cell(N,1);
data_set.measurement = cell(N,1);
data_set.t_discrete = cell(N,1);
data_set.states_discrete = cell(N,1);

for i = 1:100

% draw a sample path
system.seed = randi(intmax);
[reaction_time,state] = gillespie(system);

% generate noise observations
obs_times = [0.5*delta_t:delta_t:min(system.t_max,max(reaction_time))];
measurement = simulate_measurement(reaction_time,state,noise,obs_times,selection);

% discretize the jump process trajectories for plotting
t_discrete = [system.t_min:delta:min(system.t_max,max(reaction_time))];
[ discrete_state ] = discretize_data(reaction_time,state,t_discrete);
 
% store data
data_set.reaction_time{i} = reaction_time;
data_set.states{i} = state;
data_set.obs_times{i} = obs_times;
data_set.measurement{i} = measurement;
data_set.t_discrete{i} = t_discrete;
data_set.states_discrete{i} = discrete_state;

% plot
figure
hold on
plot(t_discrete/60,200*discrete_state(2,:),'-k','linewidth',2)
plot(t_discrete/60,discrete_state(3,:),'-b','linewidth',2)
plot(t_discrete/60,discrete_state(4,:),'-r','linewidth',2)
plot(obs_times/60,measurement,'xk','linewidth',2)
xlim([system.t_min,system.t_max/60])
ylim([0,1200])
xl = get(gca,'XAxis');
set(xl,'fontsize',14);
yl = get(gca,'YAxis');
set(yl,'fontsize',14);
xlabel('Time in min','fontsize',18)
ylabel('Abundance','fontsize',18)
title(num2str(i))
 
end

% save dataset
outfile = ['gene_expression_N_',num2str(N),'_ln.mat'];
save(outfile,'data_set')
