%% Description
% Perform exact forward backward integration for predator prey model

%% Preparations

% load required libraries
addpath('../src/ssa')
addpath('../src/models/predator_prey')
addpath('../src/models/autocatalytic_birth_death')
addpath('../src/exact_integration')

% fix random number seed for reproducibility
rng(3856820-2)


%% construct the model

% set prio rates
rates = [5e-4;1e-4;1e-4;5e-4];

Pre = [1,0; ...
    1,1; ...
    1,1; ...
    0,1];
Post = [2,0; ...
    0,1; ...
    1,2; ...
    0,0];

system = struct;
system.name = 'Predator Prey Model';
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
generator = Forward_Backward(pp_model,system);

% get initial distribution
initial_dist = generator.initial_dist(system.initial);

% solve master equation
[t,dist] = generator.integrate(initial_dist,[system.t_min,system.t_max]);

% compute marginal
population_mean = pp_model.population_mean(dist);

%% filtering and smoothing

% evaluate obs likelihood
obs_likelihood = generator.observation_likelihood(obs_model,measurement);
filter_norm = zeros(1,size(obs_likelihood,2));

% set up smoother
delta_t = 1;
smoother = Smoother(system.t_min,system.t_max,delta_t,obs_times,measurement);
smoother.set_obs_likelihood(pp_model,obs_model);

% compute filter
smoother.forward_filter(generator,initial_dist);
smoother.backward_filter(generator);

% compute the mean of filtering and smoothing distribution
filter_time = smoother.get_time();
filter_mean = pp_model.population_mean(smoother.get_forward());
smoother_mean = pp_model.population_mean(smoother.get_smoothed());

%% filtering and smoothing with autocatalytic birth death process

% reduced structure
abd_system = struct;
abd_system.rates = system.rates(1:2);
abd_system.max_state = system.max_state(1);
abd_system.initial = system.initial(1);
abd_system.t_min = 0;
abd_system.t_max = 3000;

% set up model
abd_model = Generator_ABD(abd_system.max_state);
abd_generator = Forward_Backward(abd_model,abd_system);

% get initial distribution
abd_initial_dist = abd_generator.initial_dist(abd_system.initial);

% solve master equation
[abd_t,abd_dist] = abd_generator.integrate(abd_initial_dist,[abd_system.t_min,abd_system.t_max]);

% compute marginal
abd_mean = abd_model.population_mean(abd_dist);

% set up smoother
abd_smoother = Smoother(abd_system.t_min,abd_system.t_max,delta_t,obs_times,measurement(1,:));
abd_smoother.set_obs_likelihood(abd_model,obs_model_single);

% compute filter
abd_smoother.forward_filter(abd_generator,abd_initial_dist);
abd_smoother.backward_filter(abd_generator);

% compute the mean of filtering and smoothing distribution
abd_filter_time = abd_smoother.get_time();
abd_filter_mean = abd_model.population_mean(abd_smoother.get_forward());
abd_smoother_mean = abd_model.population_mean(abd_smoother.get_smoothed());

%% plotting

log_evidence = sum(filter_norm);

% simulated trajectory and observations vs prior mean
figure
hold on
plot(t,population_mean(1,:),'--b','linewidth',2)
plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
plot(obs_times,measurement(1,:),'ok','linewidth',2)
plot(t,population_mean(2,:),'--r','linewidth',2)
plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
plot(obs_times,measurement(2,:),'xk','linewidth',2)
title('Prior Mean')
ylim([0,30])

% simulated trajectory and observations vs filter mean
figure
hold on
plot(filter_time,filter_mean(1,:),'--b','linewidth',2)
plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
plot(obs_times,measurement(1,:),'ok','linewidth',2)
plot(filter_time,filter_mean(2,:),'--r','linewidth',2)
plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
plot(obs_times,measurement(2,:),'xk','linewidth',2)
title('Filter Mean')
ylim([0,30])

% simulated trajectory and observations vs smoother mean
figure
hold on
plot(filter_time,smoother_mean(1,:),'--b','linewidth',2)
plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
plot(obs_times,measurement(1,:),'ok','linewidth',2)
plot(filter_time,smoother_mean(2,:),'--r','linewidth',2)
plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
plot(obs_times,measurement(2,:),'xk','linewidth',2)
title('Smoother Mean')
ylim([0,30])

% prey observations vs abd mean
figure
hold on
plot(abd_t,abd_mean,'--b','linewidth',2)
plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
plot(abd_filter_time,abd_filter_mean,'--r','linewidth',2)
plot(abd_filter_time,abd_smoother_mean,'--g','linewidth',2)
plot(obs_times,measurement(1,:),'ok','linewidth',2)
title('Birth Death Mean')
ylim([0,30])
