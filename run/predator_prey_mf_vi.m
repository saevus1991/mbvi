% load required libraries
addpath('../src/ssa')
addpath('../src/models/predator_prey')
addpath('../src/models/autocatalytic_birth_death')
addpath('../src/mean_field_vi')
%addpath('models/predator_prey/log_normal')

%% Description
% Apply mean-field VM to the predator prey model

%% Initialize parameters

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
generator = Forward_Backward_MF(pp_model,system);
%
%  % get initial distribution
%  initial_dist = generator.initial_dist(system.initial);
%
%  % solve master equation
%  [t,dist] = generator.integrate(initial_dist,[system.t_min,system.t_max]);
%
%  % compute marginal
%  population_mean = pp_model.population_mean(dist);

%% filtering and smoothing

% evaluate obs likelihood
obs_likelihood = generator.observation_likelihood(obs_model,measurement);
filter_norm = zeros(1,size(obs_likelihood,2));

%  % set up smoother
% delta_t = 1;
% smoother = Smoother(system.t_min,system.t_max,delta_t,obs_times,measurement);
% smoother.set_obs_likelihood(pp_model,obs_model);
%
% % compute filter
% smoother.forward_filter(generator,initial_dist);
% smoother.backward_filter(generator);
%
%  % compute the mean of filtering and smoothing distribution
%  filter_time = smoother.get_time();
%  filter_mean = pp_model.population_mean(smoother.get_forward());
%  smoother_mean = pp_model.population_mean(smoother.get_smoothed());

%% mean field variational inference

% set up the generators for both components
vi_prey = Generator_ABD(system.max_state(1));
fb_prey = Forward_Backward_MF(vi_prey,system.rates(1:2));
%generator_prey = Generator(time,generator);
vi_pred = Generator_ABD(system.max_state(2));
fb_pred = Forward_Backward_MF(vi_prey,system.rates(3:4));

% set up smoothers
delta_t = 1;
smoother_prey = Smoother_MF(system.t_min,system.t_max,delta_t,obs_times,measurement(1,:));
smoother_prey.set_obs_likelihood(vi_prey,obs_model_single);
smoother_prey.initialize_forward_generator(fb_prey);
smoother_pred = Smoother_MF(system.t_min,system.t_max,delta_t,obs_times,measurement(2,:));
smoother_pred.set_obs_likelihood(vi_pred,obs_model_single);
smoother_pred.initialize_forward_generator(fb_pred);

% initial distributions
prey_initial = fb_prey.initial_dist(system.initial(1));
pred_initial = fb_pred.initial_dist(system.initial(2));

% forward and backward for both
smoother_prey.forward_filter(fb_prey,prey_initial);
smoother_prey.backward_filter(fb_prey);
smoother_pred.forward_filter(fb_pred,pred_initial);
smoother_pred.backward_filter(fb_pred);

% get equation for plotting
filter_time = smoother_prey.get_time();
prey_filter_mean = vi_prey.population_mean(smoother_prey.get_forward());
prey_smoother_mean = vi_prey.population_mean(smoother_prey.get_smoothed());
pred_filter_mean = vi_pred.population_mean(smoother_pred.get_forward());
pred_smoother_mean = vi_pred.population_mean(smoother_pred.get_smoothed());

%% plotting

figure
hold on
plot(obs_times,measurement(1,:),'ok','linewidth',2)
plot(obs_times,measurement(2,:),'xk','linewidth',2)
plot(filter_time,prey_filter_mean,'--b','linewidth',2)
plot(filter_time,prey_smoother_mean,'-b','linewidth',2)
plot(filter_time,pred_filter_mean,'--r','linewidth',2)
plot(filter_time,pred_smoother_mean,'-r','linewidth',2)
%  log_evidence = sum(filter_norm);

%% perform vi iteration

% small constant used to replace vanishing rates in the prior process
tol = 1e-3;

% compute effective meanfield rate for prey population
mf_stats = cell(1,2);
for i = 1:size(prey_eff_rate,1)
    smoothed_dist = smoother_pred.get_smoothed(i);
    [mf_stats1,mf_stats_2] = vi_pred.mf_vi_stats(smoothed_dist,tol);
    mf_stats{1} = [ones(size(mf_stats1));mf_stats1];
    mf_stats{2} = [ones(size(mf_stats2));mf_stats2];
    smoother_prey.backward_generator{i} = fb_prey.timde_dep_generator(mf_stats);
end

% solve the backward equation
smoother_prey.filter_norm = zeros(size(smoother_prey.filter_norm));
smoother_prey.backward_filter(fb_prey);

% set up the backward generator

%
%  % simulated trajectory and observations vs prior mean
%  figure
%  hold on
%  plot(t,population_mean(1,:),'--b','linewidth',2)
%  plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
%  plot(obs_times,measurement(1,:),'ok','linewidth',2)
%  plot(t,population_mean(2,:),'--r','linewidth',2)
%  plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
%  plot(obs_times,measurement(2,:),'xk','linewidth',2)
%  title('Prior Mean')
%  ylim([0,30])
%
%   % simulated trajectory and observations vs filter mean
%  figure
%  hold on
%  plot(filter_time,filter_mean(1,:),'--b','linewidth',2)
%  plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
%  plot(obs_times,measurement(1,:),'ok','linewidth',2)
%  plot(filter_time,filter_mean(2,:),'--r','linewidth',2)
%  plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
%  plot(obs_times,measurement(2,:),'xk','linewidth',2)
%  title('Filter Mean')
%  ylim([0,30])
%
%  % simulated trajectory and observations vs smoother mean
%  figure
%  hold on
%  plot(filter_time,smoother_mean(1,:),'--b','linewidth',2)
%  plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
%  plot(obs_times,measurement(1,:),'ok','linewidth',2)
%  plot(filter_time,smoother_mean(2,:),'--r','linewidth',2)
%  plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
%  plot(obs_times,measurement(2,:),'xk','linewidth',2)
%  title('Smoother Mean')
%  ylim([0,30])
