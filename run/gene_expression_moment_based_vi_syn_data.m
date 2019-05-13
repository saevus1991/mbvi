% load required libraries
addpath('../data/synthetic_data/')
addpath('../src/moment_based_vi')
addpath('../src/models/gene_expression')
addpath('../src/models/gene_expression/central')
%addpath('models/predator_prey/log_normal')

%% Description
% Perform exact forward backward integration synthetic datasets from a
% predator prey model

%% Initialize parameters

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
 
 % prepare output data
 output = struct;
 output.evidence = cell(N,1);
 output.time = cell(N,1);
 output.control = cell(N,1);
 output.moments = cell(N,1);
 
 %% filtering and smoothing
 
 for k = 1:1
     
 % status message
 fprintf('Processing trajectory %d of %d \n',k,N)
 
 % extract data from data set
 measurement = data_set.measurement{k};
 t_discrete = data_set.t_discrete{k};
 discrete_state = data_set.states_discrete{k};
 obs_times = data_set.obs_times{k};
     


 
%% set up smoother object

% % constants
delta = 1e-4;

% set up model
model = struct;
model.num_species = 3;
model.num_reactions = 6;
model.noise = sigma;
model.initial = [1;0;0;0;0;0;0;0;0];
model.initial_rates = system.rates;
model.observed_species = [3];
model.moments = @(t,y,alpha_t,alpha) forward_equation(t,y,alpha_t,alpha);
model.constraints = @(t,y,alpha_t,alpha,moments,rates) backward_equation(t,y,alpha_t,alpha,moments,rates);
model.propensities = @(moments) propensities(moments);
model.control_gradient = @(constraints,moments) control_gradient(constraints,moments);
%additional function for the discrete version
model.moment_type = 'central';
model.operation_mode = 'smoothing';

% construct object
time_step = 0.5;
gene_expression = variational_engine(model,system.t_min,system.t_max,time_step,obs_times,measurement);


%% perform natural gradient descent

iter = 2000;
plot_iter = iter+1;
elbo_old = objective_function(gene_expression);
h = 0.01;
for i = 1:iter


[h,elbo] = natural_gradient_step(gene_expression,h,elbo_old);
%fprintf("Objective function: %f\n",elbo)
%fprintf("Step size: %f\n",h)

elbo_old = elbo;

if h < 1e-15
    disp(['Converged at iteration ',num2str(i)])
    break;
end


end

% get the moment equation
[m_t,m] = get_moments(gene_expression);
[~,u] = get_control(gene_expression);

% % mean_tr = m(:,1);
% % std_tr = sqrt(m(:,2)-m(:,1).^2);
% gene_mean_smooth = m(:,1);
% gene_err_smooth = m(:,4);
% mrna_mean_smooth =  m(:,2);
% mrna_err_smooth = sqrt(m(:,7));
% protein_mean_smooth =  m(:,3);
% protein_err_smooth = sqrt(m(:,9));

 % save stuff
 output.time{k} = m_t;
 output.evidence{k} = elbo;
 output.moments{k} = m;
 output.control{k} = u;
 
 end

 % save file
 outfile = ['../data/output/',file_path(1:end-4),'_processed.mat'];
 save(outfile,'output')