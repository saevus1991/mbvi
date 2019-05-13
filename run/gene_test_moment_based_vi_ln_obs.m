% load required libraries
addpath('../data/synthetic_data/')
addpath('../src/moment_based_vi')
addpath('../src/models/gene_expression')
addpath('../src/models/observation')

addpath('../src/models/gene_expression/central')
%addpath('models/predator_prey/log_normal')

%% Description
% Perform exact forward backward integration synthetic datasets from a
% predator prey model

%% Initialize parameters

% file with data
file_path = 'gene_expression_N_100_ln.mat';

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

%% load data

% fix cell to treat
k = 6;

% extract data from data set
measurement = data_set.measurement{k};
t_discrete = data_set.t_discrete{k};
discrete_state = data_set.states_discrete{k};
obs_times = data_set.obs_times{k};

%% prepare required structures

% set up variational model
initial = [1;0;0;0;0;0;0;0;0];
initial_rates = system.rates;
tspan = [system.t_min,system.t_max];
model = GeneExpressionCentral(initial,initial_rates,tspan);

% set up observation model
num_species = 3;
observed_species = 3;
obs_model = LognormObs(sigma,num_species,observed_species);

%% 

% set up options
options = struct;
options.operation_mode = 'smoothing';

% set up vi engine
time_step = 0.5;
gene_expression = variational_engine(model,obs_model,time_step,obs_times,measurement,options);

% %% set up smoother object
% 
% % % constants
% delta = 1e-4;
% 
% % set up model
% model = struct;
% model.num_species = 3;
% model.num_reactions = 6;
% model.noise = sigma;
% model.initial = [1;0;0;0;0;0;0;0;0];
% model.initial_rates = system.rates;
% model.observed_species = [3];
% model.moments = @(t,y,alpha_t,alpha) forward_equation(t,y,alpha_t,alpha);
% model.constraints = @(t,y,alpha_t,alpha,moments,rates) backward_equation(t,y,alpha_t,alpha,moments,rates);
% model.propensities = @(moments) propensities(moments);
% model.control_gradient = @(constraints,moments) control_gradient(constraints,moments);
% %additional function for the discrete version
% model.moment_type = 'central';
% model.operation_mode = 'smoothing';
% 
% % construct object
% time_step = 0.5;
% gene_expression = variational_engine(model,system.t_min,system.t_max,time_step,obs_times,measurement);
% 

%% perform natural gradient descent

iter = 200;
plot_iter = 10;
elbo_old = objective_function(gene_expression);
h = 0.001;
obs_model.reset_type = 'full';
for i = 1:iter
    
    
    [h,elbo] = gradient_step(gene_expression,h,elbo_old);
    fprintf("Objective function: %f\n",elbo)
    %fprintf("Step size: %f\n",h)
    
    elbo_old = elbo;
    
    if h < 1e-15
        disp(['Converged at iteration ',num2str(i)])
        break;
    end
    
    %change from reduced to full gradient
    if i == 10
        obs_model.reset_type = 'full';
    end
  
    % get the moment equation
[m_t,m] = get_moments(gene_expression);
[~,u] = get_control(gene_expression);

% mean_tr = m(:,1);
% std_tr = sqrt(m(:,2)-m(:,1).^2);
gene_mean_smooth = m(:,1);
gene_err_smooth = sqrt(m(:,4));
mrna_mean_smooth =  m(:,2);
mrna_err_smooth = sqrt(m(:,7));
protein_mean_smooth =  m(:,3);
protein_err_smooth = sqrt(m(:,9));

if rem(i,plot_iter) == 0
figure
hold on
plot(t_discrete,200*discrete_state(2,:),':k','linewidth',1)
p1 = plot(m_t,200*gene_mean_smooth,'-k','linewidth',2);
plot(t_discrete,discrete_state(3,:),':b','linewidth',1)
p2 = plot(m_t,mrna_mean_smooth,'-b','linewidth',2);
plot(m_t,mrna_mean_smooth+mrna_err_smooth,'-b','linewidth',1)
plot(m_t,mrna_mean_smooth-mrna_err_smooth,'-b','linewidth',1)
plot(t_discrete,discrete_state(4,:),':r','linewidth',1)
p3 = plot(m_t,protein_mean_smooth,'-r','linewidth',2);
plot(m_t,protein_mean_smooth+protein_err_smooth,'-r','linewidth',1)
plot(m_t,protein_mean_smooth-protein_err_smooth,'-r','linewidth',1)
p4 = plot(obs_times,measurement,'xk','linewidth',2);
xtick = get(gca,'XTickLabel');
set(gca,'linewidth',2);
set(gca,'XTickLabel',xtick,'fontsize',14) % FontName','Times'
xlabel('Time','interpreter','latex','fontsize',18')
ylabel('Abundance','interpreter','latex','fontsize',18')
legend([p1,p2,p3,p4],{'Gene','mRNA','Protein','Observation'},'fontsize',18,'interpreter','latex')
xlim([min(m_t),max(m_t)])
ylim([0,600])
end

% figure
% for j = 1:6
% subplot(3,2,j)
% plot(m_t,u(:,j)/initial_rates(j))
% end

    
end



% 
% % save file
% outfile = ['../data/output/',file_path(1:end-4),'_processed.mat'];
% save(outfile,'output')