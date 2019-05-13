%% Description
% Perform exact forward backward integration synthetic datasets from a
% predator prey model (does not work because of to large state space)


%% Preparations

% load required libraries
addpath('../data/synthetic_data/')
addpath('../src/exact_integration')
addpath('../src/models/gene_expression')
%addpath('models/predator_prey/log_normal')

% file with data
file_path = 'gene_expression_N_100.mat';

% load the dataset
load(file_path);

% number of files
N = length(data_set.measurement);

%% load the model

system = data_set.system;

% Gaussian observation noise
 sigma = data_set.sigma;
 noise = data_set.noise_model;
 
 obs_model = @(x,y) -log(2*pi*sigma^2)-sum(log(y))-0.5*sum((log(y)-log(x+1)).^2)/sigma^2;

% constants
delta = 1;
delta_t = 300;

 %% integration with state space for prior
 
 % set up integrator object
 max_state = [1;200;800];
 GE_model = Generator_GE(max_state);
 generator = Forward_Backward(GE_model,system);
 
 % get initial distribution
 initial_dist = generator.initial_dist(system.initial);
 
 % solve master equation
 [t,dist] = generator.integrate(initial_dist,[system.t_min,system.t_max]);
 
 % compute marginal
 population_mean = pp_model.population_mean(dist);
 
 % prepare output data
 output = struct;
 output.evidence = cell(N,1);
 output.filter_time = cell(N,1);
 output.filter_stats = cell(N,2);
 output.smoother_stats = cell(N,2);
 
 %% filtering and smoothing
 
 for i = 1:0
     
 % status message
 fprintf('Processing trajectory %d of %d \n',i,N)
 
 % extract data from data set
 measurement = data_set.measurement{i};
 t_discrete = data_set.t_discrete{i};
 discrete_state = data_set.states_discrete{i};
 obs_times = data_set.obs_times{i};
     
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
 filter_std = pp_model.population_std(smoother.get_forward(),filter_mean);
 smoother_mean = pp_model.population_mean(smoother.get_smoothed());
 smoother_std = pp_model.population_std(smoother.get_smoothed(),smoother_mean);

 % save stuff
 output.filter_time{i} = filter_time;
 output.evidence{i} = sum(smoother.filter_norm);
 output.filter_stats{i,1} = filter_mean;
 output.filter_stats{i,2} = filter_std;
 output.smoother_stats{i,1} = smoother_mean;
 output.smoother_stats{i,2} = smoother_std;
 
%  % simulated trajectory and observations vs smoother mean
%  figure
%  hold on
%  plot(filter_time,smoother_mean(1,:),'--b','linewidth',2)
%  plot(filter_time,smoother_mean(1,:)+smoother_std(1,:),':b','linewidth',1)
%  plot(filter_time,smoother_mean(1,:)-smoother_std(1,:),':b','linewidth',1)
%  plot(t_discrete,discrete_state(1,:),'-b','linewidth',2)
%  plot(obs_times,measurement(1,:),'ok','linewidth',2)
%  plot(filter_time,smoother_mean(2,:),'--r','linewidth',2)
%  plot(filter_time,smoother_mean(2,:)+smoother_std(2,:),':r','linewidth',1)
%  plot(filter_time,smoother_mean(2,:)-smoother_std(2,:),':r','linewidth',1)
%  plot(t_discrete,discrete_state(2,:),'-r','linewidth',2)
%  plot(obs_times,measurement(2,:),'xk','linewidth',2)
%  title('Smoother Mean')
%  ylim([0,30])
 
 end

 % save file
 outfile = ['../data/output/',file_path(1:end-4),'_processed.mat'];
 save(outfile,'output')