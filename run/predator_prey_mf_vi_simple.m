%% Description
% Apply mean-field VM to the predator prey model
% do not use full sparse matrix objective but direct implementation via
% birth death processes

%% Preparations

% load required libraries
addpath('../src/ssa')
addpath('../src/models/predator_prey')
addpath('../src/models/autocatalytic_birth_death')
addpath('../src/mean_field_vi')
%addpath('models/predator_prey/log_normal')

% fix random number seed for reproducibility
rng(3856820-2)

% construct priors
rates = [5e-4;1e-4;1e-4;5e-4];


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
obs_model_single = @(x,y) -0.5*log(2*pi*sigma^2)-log(y)-0.5*(log(y)-log(x+1)).^2/sigma^2;

% constants
delta = 1;
delta_t = 300;

%% simulation

% observed species
selection = [1;2];

% draw a sample path
[reaction_time,state] = gillespie(system);

% generate noise observations
obs_times = [0.5*delta_t:delta_t:min(system.t_max,max(reaction_time))];
measurement = simulate_measurement(reaction_time,state,noise,obs_times,selection);

% discretize the jumo process trajectories for plotting
t_discrete = [system.t_min:delta:min(system.t_max,max(reaction_time))];
[ discrete_state ] = discretize_data(reaction_time,state,t_discrete);

rate_equation = @(t,m) [rates(1)*m(1)-rates(2)*m(1)*m(2);rates(3)*m(1)*m(2)-rates(4)*m(2)];
[rre_t,rre_state] = ode45(rate_equation,t_discrete,system.initial);


%% integration with state space for prior

% set up integrator object
pp_model = Generator_PP(system.max_state);


%% mean field variational inference

% set up the generators for both components
vi_prey = Generator_ABD(system.max_state(1));
vi_pred = Generator_ABD(system.max_state(2));

% solve components forward in time
vi_prey = 



%% plotting

figure
hold on
plot(obs_times,measurement(1,:),'ok','linewidth',2)
plot(obs_times,measurement(2,:),'xk','linewidth',2)
plot(t_discrete,discrete_state(1,:),'--b','linewidth',2)
plot(t_discrete,discrete_state(2,:),'--r','linewidth',2)
% plot(filter_time,prey_filter_mean,'--b','linewidth',2)
% plot(filter_time,prey_smoother_mean,'-b','linewidth',2)
% plot(filter_time,pred_filter_mean,'--r','linewidth',2)
% plot(filter_time,pred_smoother_mean,'-r','linewidth',2)


