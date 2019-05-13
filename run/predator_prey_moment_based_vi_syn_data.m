% load required libraries
addpath('../data/synthetic_data/')
addpath('../src/models/predator_prey')
addpath('../src/models/predator_prey/log_normal')
addpath('../src/moment_based_vi')
%addpath('models/predator_prey/log_normal')

%% Description
% Perform exact forward backward integration synthetic datasets from a
% predator prey model

%% Initialize parameters

% file with data
file_path = 'predator_prey_N_100.mat';

% load the dataset
load(file_path);

% number of files
N = length(data_set.measurement);

%% load the model

system = data_set.system;

% Gaussian observation noise
sigma = 3;%data_set.sigma;
noise = data_set.noise_model;
obs_model = @(x,y) -log(2*pi*sigma^2)-sum(log(y))-0.5*sum((log(y)-log(x+1)).^2)/sigma^2;

% observed speciees
selection = [1;2];

% constants
delta = 1;
delta_t = 300;

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

% prepare output data
output = struct;
output.evidence = cell(N,1);
output.filter_time = cell(N,1);
output.filter_stats = cell(N,2);
output.smoother_stats = cell(N,2);

%% filtering and smoothing

for k = 5:5
    
    % status message
    fprintf('Processing trajectory %d of %d \n',k,N)
    
    % extract data from data set
    measurement = data_set.measurement{k};
    t_discrete = data_set.t_discrete{k};
    discrete_state = data_set.states_discrete{k};
    obs_times = data_set.obs_times{k};
    
    
    %% set up smoother object
    
    % set up model
    model = struct;
    model.num_species = 2;
    model.num_reactions = 4;
    model.noise = sigma;
    model.initial = [20;5;0;0;0];
    model.initial_rates = system.rates;
    model.observed_species = selection;
    model.moments = @(t,y,alpha_t,alpha) forward_equation(t,y,alpha_t,alpha);
    model.constraints = @(t,y,alpha_t,alpha,moments,rates) backward_equation(t,y,alpha_t,alpha,moments,rates);
    model.propensities = @(moments) propensities(moments);
    model.control_gradient = @(constraints,moments) control_gradient(constraints,moments);
    %additional function for the discrete version
    model.moment_type = 'central';
    model.operation_mode = 'smoothing';
    model.backward_jump = @(m,y) backward_jump_ln3(m,y,0.15);
    %model.backward_jump = @(m,y) backward_jump_gauss(m,y,3);
    
    % construct object
    time_step = 0.5;
    gene_expression = variational_engine(model,system.t_min,system.t_max,time_step,obs_times,measurement);
    
    % get the moment equation
    [m_t,m] = get_moments(gene_expression);
    
    % mean_tr = m(:,1);
    % std_tr = sqrt(m(:,2)-m(:,1).^2);
    prey_mean = m(:,1);
    prey_err = sqrt(m(:,3));
    predator_mean =  m(:,2);
    predator_err = sqrt(m(:,5));
    
    
    figure
    hold on
    plot(m_t/60,prey_mean,'-b','linewidth',2)
    plot(m_t/60,predator_mean,'-r','linewidth',2)
    plot(t_discrete/60,discrete_state(1,:),'--b','linewidth',1)
    plot(m_t/60,prey_mean+prey_err,':b','linewidth',1)
    plot(t_discrete/60,discrete_state(2,:),'--r','linewidth',1)
    plot(m_t/60,predator_mean+predator_err,':r','linewidth',1)
    plot(obs_times/60,measurement(1,:),'xk','linewidth',2)
    if size(measurement,1) == 2
        plot(obs_times/60,measurement(2,:),'ok','linewidth',2)
    end
    plot(m_t/60,prey_mean-prey_err,':b','linewidth',1)
    plot(m_t/60,predator_mean-predator_err,':r','linewidth',1)
    legend({'Prey','Predator'},'fontsize',18)
    xlim([system.t_min,system.t_max/60])
    %ylim([0,50])
    xl = get(gca,'XAxis');
    set(xl,'fontsize',14);
    yl = get(gca,'YAxis');
    set(yl,'fontsize',14);
    xlabel('Time in min','fontsize',18)
    ylabel('Abundance','fontsize',18)
    %print('extended_gene_expression_smoothing_prior','-depsc')
    
    
    %% rund gradient based smoothing
    
    iter = 1000;
    elbo_old = objective_function(gene_expression);
    h = 0.01;
    for i = 1:iter
        
        %elbo = partial_projected_gradient_descent(gene_expression,[1;2])
        %elbo = projected_gradient_descent(gene_expression)
        %elbo = exp_gradient_descent(gene_expression)
        [h,elbo] = natural_gradient_step(gene_expression,h,elbo_old);
        fprintf("Objective function: %f\n",elbo)
        fprintf("Step size: %f\n",h)
        
        elbo_old = elbo;
        
        if h < 1e-15
            disp(['Converged at iteration ',num2str(i)])
            break;
        end
        
        
    end
    
    % get the moment equation
    [m_t,m] = get_moments(gene_expression);
    
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
    
    %  % save stuff
    %  output.filter_time{k} = filter_time;
    %  output.evidence{k} = sum(smoother.filter_norm);
    %  output.filter_stats{k,1} = filter_mean;
    %  output.filter_stats{k,2} = filter_std;
    %  output.smoother_stats{k,1} = smoother_mean;
    %  output.smoother_stats{k,2} = smoother_std;
    
end

% save file
%outfile = ['../data/output/',file_path(1:end-4),'_processed.mat'];
%save(outfile,'output')



function jump = backward_jump_gauss(m,y,sigma)

jump = zeros(size(m,1),size(y,2));
jump(1,:) = (y(1,:)-m(1,:))/sigma^2;
jump(2,:) = (y(2,:)-m(2,:))/sigma^2;
jump(3,:) = -0.5/sigma^2;
jump(5,:) = -0.5/sigma^2;

end

function jump = backward_jump_ln(m,y,sigma)

m1_eff = 1+m(1,:);
m2_eff = 1+m(2,:);
jump = zeros(size(m,1),size(y,2));
jump(1,:) = -1./m1_eff-2*m1_eff./(m(3,:)-m1_eff.^2)+2*(log(m1_eff)-0.5*log(m(3,:)))./(m1_eff)-2*log(y(1,:))./m1_eff;
jump(2,:) = -1./m2_eff-2*m2_eff./(m(5,:)-m2_eff.^2)+2*(log(m2_eff)-0.5*log(m(5,:)))./(m2_eff)-2*log(y(2,:))./m2_eff;
jump(3,:) = 1./(m(3,:)-m1_eff.^2)-(log(m1_eff)-0.5*log(m(3,:)))./m(3,:)+log(y(1,:))./m(3,:);
jump(5,:) = 1./(m(5,:)-m2_eff.^2)-(log(m2_eff)-0.5*log(m(5,:)))./m(5,:)+log(y(1,:))./m(5,:);

jump = -0.5*jump/sigma^2;

end

function jump = backward_jump_ln2(m,y,sigma)

m1_eff = 1+m(1,:);
m2_eff = 1+m(2,:);
jump = zeros(size(m,1),size(y,2));
jump(1,:) = -2./m1_eff+2*m1_eff./(m(3,:)+m1_eff.^2)+2*(2./m1_eff-m1_eff./(m1_eff.^2+m(3,:))).*(2*log(m1_eff)-0.5*log(m1_eff.^2+m(3,:)))-4*(1./m1_eff-0.5*m1_eff./(m1_eff.^2+m(3,:))).*log(y(1,:));
jump(2,:) = -2./m2_eff+2*m2_eff./(m(5,:)+m2_eff.^2)+2*(2./m2_eff-m2_eff./(m2_eff.^2+m(5,:))).*(2*log(m2_eff)-0.5*log(m2_eff.^2+m(5,:)))-4*(1./m2_eff-0.5*m2_eff./(m2_eff.^2+m(5,:))).*log(y(2,:));
jump(3,:) = 1./(m(3,:)+m1_eff.^2)-(2*log(m1_eff)-0.5*log(m(3,:)+m1_eff.^2))./(m(3,:)+m1_eff.^2)+log(y(1,:))./(m(3,:)+m1_eff.^2);
jump(5,:) = 1./(m(5,:)+m2_eff.^2)-(2*log(m2_eff)-0.5*log(m(5,:)+m2_eff.^2))./(m(5,:)+m2_eff.^2)+log(y(2,:))./(m(5,:)+m2_eff.^2);

jump = -0.5*jump/sigma^2;

end

function jump = backward_jump_ln3(m,y,sigma)

m1_eff = 1+m(1,:);
m2_eff = 1+m(2,:);
m3_eff = m1_eff.^2+m(3,:);
m5_eff = m2_eff.^2+m(5,:);
jump = zeros(size(m,1),size(y,2));
jump(1,:) = -2./m1_eff+2*m1_eff./m3_eff+2*(2./m1_eff-m1_eff./m3_eff).*(2*log(m1_eff)-0.5*log(m3_eff))-4*(1./m1_eff-0.5*m1_eff./m3_eff).*log(y(1,:));
jump(2,:) = -2./m2_eff+2*m2_eff./m5_eff+2*(2./m2_eff-m2_eff./m5_eff).*(2*log(m2_eff)-0.5*log(m5_eff))-4*(1./m2_eff-0.5*m2_eff./m5_eff).*log(y(2,:));
jump(3,:) = (1-(2*log(m1_eff)-0.5*log(m3_eff))+log(y(1,:)))./m3_eff;
jump(5,:) = (1-(2*log(m2_eff)-0.5*log(m5_eff))+log(y(2,:)))./m5_eff;

jump = -0.5*jump/sigma^2;

end