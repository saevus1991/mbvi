% load required libraries
addpath('../data/synthetic_data/')
addpath('../src/models/predator_prey/log_normal')
addpath('../src/models/observation')
addpath('../src/moment_based_vi')
%addpath('models/predator_prey/log_normal')

%% Description
% Perform exact forward backward integration synthetic datasets from a
% predator prey model

%% Initialize parameters

% file with data
file_path = 'predator_prey_gauss_N_100.mat';

% load the dataset
load(file_path);

% number of files
N = length(data_set.measurement);

%% Preparations

% load system
system = data_set.system;

% observed speciees
selection = [1;2];

% constants
delta = 1;
delta_t = 300;

% set up variational model
initial = [20;5;0;0;0];
initial_rates = system.rates;
tspan = [system.t_min,system.t_max];
model = PredatorPreyLognormal(initial,initial_rates,tspan);

% set up observation model
sigma = [data_set.sigma;data_set.sigma];
num_species = 2;
observed_species = [1;2];
obs_model = MultiGaussObs(sigma,num_species,observed_species);

% set up options
options = struct;
options.operation_mode = 'smoothing';

% prepare output data
output = struct;
output.evidence = cell(N,1);
output.time = cell(N,1);
output.control = cell(N,1);
output.moments = cell(N,1);

%% filtering and smoothing

for k = 23:23
    
    % status message
    fprintf('Processing trajectory %d of %d \n',k,N)
    
    % extract data from data set
    measurement = data_set.measurement{k};
    t_discrete = data_set.t_discrete{k};
    discrete_state = data_set.states_discrete{k};
    obs_times = data_set.obs_times{k};
    
    % set up vi engine
    time_step = 0.5;
    predator_prey = variational_engine(model,obs_model,time_step,obs_times,measurement,options);  
    
    %% rund gradient based smoothing
    
    iter = 2000;
    elbo_old = objective_function(predator_prey);
    h = 0.01;
    for i = 1:iter
        
        %elbo = partial_projected_gradient_descent(predator_prey,[1;2])
        %elbo = projected_gradient_descent(predator_prey)
        %elbo = exp_gradient_descent(predator_prey)
        [h,elbo] = natural_gradient_step(predator_prey,h,elbo_old);
        %         fprintf("Objective function: %f\n",elbo)
        %         fprintf("Step size: %f\n",h)
        
        elbo_old = elbo;
        
        if h < 1e-15
            disp(['Converged at iteration ',num2str(i)])
            break;
        end
        
        
    end
    
    % get the moment equation
    [m_t,m] = get_moments(predator_prey);
    [~,u] = get_control(predator_prey);
    
    prey_mean_smooth = m(:,1);
    prey_err_smooth = sqrt(m(:,3));
    predator_mean_smooth =  m(:,2);
    predator_err_smooth = sqrt(m(:,5));
    
    figure
    hold on
    plot(m_t/60,prey_mean_smooth,'-b','linewidth',2)
    plot(m_t/60,predator_mean_smooth,'-r','linewidth',2)
    plot(t_discrete/60,discrete_state(1,:),'--b','linewidth',1)
    plot(m_t/60,prey_mean_smooth+prey_err_smooth,':b','linewidth',1)
    plot(t_discrete/60,discrete_state(2,:),'--r','linewidth',1)
    plot(m_t/60,predator_mean_smooth+predator_err_smooth,':r','linewidth',1)
    plot(obs_times/60,measurement(1,:),'xk','linewidth',2)
    if size(measurement,1) == 2
        plot(obs_times/60,measurement(2,:),'ok','linewidth',2)
    end
    plot(m_t/60,prey_mean_smooth-prey_err_smooth,':b','linewidth',1)
    plot(m_t/60,predator_mean_smooth-predator_err_smooth,':r','linewidth',1)
    legend({'Prey','Predator'},'fontsize',18)
    xlim([system.t_min,system.t_max/60])
    %ylim([0,50])
    xl = get(gca,'XAxis');
    set(xl,'fontsize',14);
    yl = get(gca,'YAxis');
    set(yl,'fontsize',14);
    xlabel('Time in min','fontsize',18)
    ylabel('Abundance','fontsize',18)
    
    % save stuff
    output.time{k} = m_t;
    output.evidence{k} = elbo;
    output.moments{k} = m;
    output.control{k} = u;
    
end

% save file
% outfile = ['../data/output/',file_path(1:end-4),'_processed.mat'];
% save(outfile,'output')


