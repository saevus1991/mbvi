% load required libraries
addpath('../data/synthetic_data/')
addpath('../data/output')
addpath('../src/models/gene_expression')
addpath('../src/models/gene_expression/central')
%addpath('models/predator_prey/log_normal')

%% Description
% Produce plot from a gene expression model with the moment based smoothing
% predator prey model

%% Load dataset

% load the simulations dataset
file_path = 'gene_expression_N_100.mat';
sim_data = load(file_path);
sim_data = sim_data.data_set;

% load the vi dataset
file_path = 'gene_expression_N_100_processed.mat';
mmvi_data = load(file_path);
mmvi_data = mmvi_data.output;

% number of files
N = length(sim_data.measurement);

%% load the model

system = sim_data.system;

% Gaussian observation noise
sigma = sim_data.sigma;
noise = sim_data.noise_model;

% pick file
k = 5;

% extract data from data set
measurement = sim_data.measurement{k};
t_discrete = sim_data.t_discrete{k};
discrete_state = sim_data.states_discrete{k};
obs_times = sim_data.obs_times{k};

% convert the moments
m = mmvi_data.moments{k};
m_t = mmvi_data.time{k};
gene_mean_smooth = m(:,1);
gene_err_smooth = sqrt(m(:,4));
mrna_mean_smooth =  m(:,2);
mrna_err_smooth = sqrt(m(:,7));
protein_mean_smooth =  m(:,3);
protein_err_smooth = sqrt(m(:,9));

%% plotting

%  % simulated mean vs rre mean
%  figure
%  hold on
%  plot(t_discrete,moments(2,:),'-b','linewidth',2)
%  plot(rre_t,rre_state(:,2),'--b','linewidth',2)
%  plot(t_discrete,moments(3,:),'-r','linewidth',2)
%  plot(rre_t,rre_state(:,3),'--r','linewidth',2)
%  %ylim([0,50])

% simulated mean vs truncated mean
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

% % reduce whitespace
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

% save
filename = 'gene_expression_moment_based_vi';
print(filename,'-dpdf')