% load required libraries
addpath('../../src/ssa')
addpath('../../src/models/predator_prey')
%addpath('models/predator_prey/log_normal')

%% Description
% Create synthetic dataset for predator prey model

%% Initialize parameters

% fix random number seed for reproducibility
rng(3856820-2)

% construct priors
rates = [5e-4;1e-4;1e-4;5e-4];

% number of trajectories
N = 100;


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
 sigma = 0.15;
 noise = @(x) lognrnd(log(x+1),sigma);
 
 obs_model = @(x,y) -log(2*pi*sigma^2)-sum(log(y))-0.5*sum((log(y)-log(x+1)).^2)/sigma^2;

% constants
delta = 1;
delta_t = 300;

%% simulation

% observed species
selection = [1;2];

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

for i = 1:N

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
 
end

% save dataset
outfile = ['predator_prey_N_',num2str(N),'.mat'];
save(outfile,'data_set')
