%% Description

% run birth death process VI with a simple Gaussian noise endpoint
% observation

%% Preparations

%rng(4856823)

% load required libraries
addpath('../src/moment_based_vi')
addpath('../src/models/birth_death/central')
addpath('../src/models/observation')

% simulation parameters
T_max = 60;            % runtime of simulation
d1 = 5;                % birth rate
d2 = 0.1;                % death rate per individual

%% construct the model

% define a the birth death process model
Pre = [0; ...
      1];
Post = [1; ...
       0];
rates = [d1;d2];

a = [1.5;1.5];
b = a./rates;

system = struct;
system.name = 'Birth Death Process';
system.Pre = Pre;
system.Post = Post;
system.S = Post-Pre;
system.rates = rates;
system.initial = [0];
system.t_min = 0;
system.t_max = T_max;
system.num_species = size(system.Pre,2);
system.num_reactions = size(system.Pre,1);
system.seed = randi(intmax);

% Gaussian observation noise
sigma = 15;
noise = @(x) normrnd_2(x,sigma);

%% aritificial data
sample_time = 50;
measurement = 0;


%% compute exact solution

% compute the prior solution
alpha_t = linspace(0,50,1000)';
alpha = zeros(length(alpha_t),2);
alpha(:,1) = d1;
alpha(:,2) = d2;

forward_update = @(t,y) forward_equation(t,y,alpha_t,alpha);
[~,exact_moment] = ode15s(forward_update,alpha_t',[0;0]);

prior_mean = exact_moment(:,1);
prior_std = sqrt(exact_moment(:,2));

%prior_mean = [prior_mean;zeros(2,1)];
%prior_std = [prior_std;zeros(2,1)];

% set up the constraints
alpha_t = linspace(0,50,1000)';
alpha = zeros(length(alpha_t),2);
alpha(:,1) = d1*(1-exp(-d2*(50-alpha_t)));
alpha(:,2) = d2./(1-exp(-d2*(50-alpha_t)));

forward_update = @(t,y) forward_equation(t,y,alpha_t,alpha);
[~,exact_moment] = ode15s(forward_update,alpha_t',[0;0]);

exact_mean = exact_moment(:,1);
exact_std = sqrt(exact_moment(:,2));

exact_mean = [exact_mean;zeros(2,1)];
exact_std = [exact_std;zeros(2,1)];


%% set up smoother object

% % constants
delta = 1e-4;

% set up model
% model = struct;
% model.num_species = 1;
% model.num_reactions = 2;
% model.alpha = a;
% model.beta = b;
% model.noise = sigma;
% model.observed_species = [1];
% model.moments = @(t,y,alpha_t,alpha) forward_equation(t,y,alpha_t,alpha);
% model.constraints = @(t,y,alpha_t,alpha,rates) backward_equation(t,y,alpha_t,alpha,rates);
% model.propensities = @(moments) propensities(moments);
% model.control_gradient = @(constraints,moments) control_gradient(constraints,moments);
% model.moment_type = 'central';

% construct gene expression model
initial = [0;0];
initial_rates = rates;
tspan = [0 T_max];
model = BirthDeathCentral(initial,initial_rates,tspan);

% set up observation model
num_species = 1;
observed_species = 1;
obs_model = GaussObs(sigma,num_species,observed_species);

% set up options
options = struct;
options.operation_mode = 'smoothing';
options.control_bounds = [1e-6,1e6];


% construct object
time_step = 0.1;
birth_death = variational_engine(model,obs_model,time_step,sample_time,measurement,options);




%% test the forward backward sweep method for central moments


obs_model.sigma = 10;

iter = 100;
elbo_old = Inf;
h = 0.01;

for i = 1:iter
    % natural gradient step
    [h,elbo] = gradient_step(birth_death,h,elbo_old);
    fprintf("Objective function: %f\n",elbo)
    %fprintf("Step size: %f\n",h)
    
    elbo_old = elbo;
    
    if h < 1e-15
        disp(['Converged at iteration ',num2str(i)])
        break;
    end
end

% get the moment equation
[m_t,m] = get_moments(birth_death);

noise1_t = m_t;
noise1_mean =  m(:,1);
noise1_std = sqrt(m(:,2));


obs_model.sigma = 5;
elbo_old = inf;
h = 0.01;

for i = 1:iter
    % natural gradient step
    [h,elbo] = gradient_step(birth_death,h,elbo_old);
    fprintf("Objective function: %f\n",elbo)
    %fprintf("Step size: %f\n",h)
    
    elbo_old = elbo;
    
    if h < 1e-15
        disp(['Converged at iteration ',num2str(i)])
        break;
    end
end

% get the moment equation
[m_t,m] = get_moments(birth_death);

noise2_t = m_t;
noise2_mean =  m(:,1);
noise2_std = sqrt(m(:,2));


obs_model.sigma = 2;
elbo_old = inf;
h = 0.01;

for i = 1:iter
    % natural gradient step
    [h,elbo] = gradient_step(birth_death,h,elbo_old);
    fprintf("Objective function: %f\n",elbo)
    %fprintf("Step size: %f\n",h)
    
    elbo_old = elbo;
    
    if h < 1e-15
        disp(['Converged at iteration ',num2str(i)])
        break;
    end
end

noise3_t = m_t;
noise3_mean =  m(:,1);
noise3_std = sqrt(m(:,2));


figure
hold on
plot(alpha_t,prior_mean,'-b','linewidth',2)
plot(noise1_t,noise1_mean,'-g','linewidth',2)
plot(noise2_t,noise2_mean,'-m','linewidth',2)
plot(alpha_t,exact_mean,'-r','linewidth',2)
plot(alpha_t,exact_mean+exact_std,':r','linewidth',1)
plot(alpha_t,exact_mean-exact_std,':r','linewidth',1)
plot(alpha_t,prior_mean+prior_std,':b','linewidth',1)
plot(alpha_t,prior_mean-prior_std,':b','linewidth',1)
plot(noise1_t,noise1_mean+noise1_std,':g','linewidth',1)
plot(noise1_t,noise1_mean-noise1_std,':g','linewidth',1)
plot(noise2_t,noise2_mean+noise2_std,':m','linewidth',1)
plot(noise2_t,noise2_mean-noise2_std,':m','linewidth',1)
legend({'Prior',' $\sigma=10$','$\sigma=5$','Noiseless'},'fontsize',18,'interpreter','latex','location','south')
xlim([system.t_min,50])
ylim([0,60])
xtick = get(gca,'XTickLabel');
set(gca,'linewidth',2);
set(gca,'XTickLabel',xtick,'fontsize',14) 
xlabel('Time in s','fontsize',18,'interpreter','latex')
ylabel('Abundance','fontsize',18,'interpreter','latex')


print('birth_death_endpoint','-dpdf')



