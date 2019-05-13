%% Description
% test moment based vi with product closure for the SIS model

%% Preparations

rng(2389514379)

addpath('../src/models/agent_based/sis')
addpath('../src/models/observation')
addpath('../src/ssa')
addpath('../src/moment_based_vi')

% specify graph by adjecency matrix
A = [0,1,0,0,0,0;
    1,0,1,1,0,0;
    0,1,0,0,0,0;
    0,1,0,0,1,1;
    0,0,0,1,0,0;
    0,0,0,1,0,0];

% initial state
x = [1;0;0;0;0;0];
rates = [0.04;0.0];

% set up a system structure for stochastic simulation
system.initial = x;
system.rates = rates;      % infection rate, recovery rate
system.graph = A;
system.dynamcis = 'SIS';
system.t_min = 0;
system.t_max = 100;

% observations
obs_times = [25,75];
sigma = 0.05;
selection = 1:6;
noise = @(x) normrnd(x,sigma);

% set up integrator object
integrator = Generator_SIS(system);

% set up variational model
initial = [1-x,x]';
initial = initial(:);
initial(2:2:end) = x;
initial_rates = repmat([0;rates;0],[length(x),1]);
model = SIS_Model(initial,initial_rates,A,[system.t_min,system.t_max],2);

% set up observation model
num_species = 3;
observed_species = 3;
obs_model = GaussObsBernoulli(sigma,model.num_nodes);

% set up options
options = struct;
options.operation_mode = 'smoothing';

% simulate synthetic data
t_sample = linspace(system.t_min,system.t_max,1000);
[times,states] = ssa_sis(system);
discrete_data = discretize_data([system.t_min,times],[system.initial,states],t_sample);
measurement = simulate_measurement(times,states,noise,obs_times,selection);

% preapre vi engine
time_step = 0.5;
vi_engine = variational_engine(model,obs_model,time_step,obs_times,measurement,options);

%% perform computations

% compute the forward mean field equation using the sis model
alpha_t = [system.t_min;system.t_max];
alpha = [0;rates;0];
alpha = repmat(alpha,[model.num_nodes,1]);
alpha = [alpha,alpha]';
fun = @(t,m)model.moments(t,m,alpha_t,alpha);
[t_ode,states_ode] = ode45(fun,[system.t_min;system.t_max-1e-3],initial);
[t_vi,states_vi] = get_moments(vi_engine);

% get initial distribution
initial_dist = integrator.initial_dist(system.initial);

% solve master equation
[t,dist] = integrator.integrate(initial_dist,system.t_min,system.t_max);

% compute marginal
marginal_dist = integrator.marginalize(dist);

% compute meanfield solutions of the marginals
ode_fun = @(t,m) simple_moment_equation(t,m,system);
[t_mf,marginal_mf] = ode45(ode_fun,[system.t_min,system.t_max],system.initial);
marginal_mf = marginal_mf';

% compute a sample average using simulations
% N = 1000;
% marginal_sim = zeros(length(system.initial),N);
% for i = 1:N
%     [times,states] = ssa_sis(system);
%     discrete_data = discretize_data([system.t_min,times],[system.initial,states],t_sample);
%     marginal_sim = marginal_sim+discrete_data;
% end
% marginal_sim = marginal_sim/N;

%% perform natural gradient descent

iter = 5;
plot_iter = iter+1;
elbo_old = objective_function(vi_engine);
h = 0.01;
for j = 1:iter
    
    
    [h,elbo] = gradient_step(vi_engine,h,elbo_old);
    fprintf("Objective function: %f\n",elbo)
    %fprintf("Step size: %f\n",h)
    
    elbo_old = elbo;
    
    if h < 1e-15
        disp(['Converged at iteration ',num2str(i)])
        break;
    end

% get data
[t_vi,states_vi] = get_moments(vi_engine);

% plot marginals
figure
for i = 1:6
    subplot(3,2,i)
    hold on
    plot(t,marginal_dist(i,:),'-k','linewidth',2)
    plot(t_mf,marginal_mf(i,:),'-g','linewidth',2)
    plot(t_sample,discrete_data(i,:),'--r','linewidth',2)
    plot(obs_times,measurement(i,:),'xk','linewidth',2)
    plot(t_vi,states_vi(:,2*i),'-b','linewidth',2)
    ylabel(num2str(i),'fontsize',14)
    xlabel('Time','fontsize',14)
    ylim([-0.1,1.1])
end


    
    
end



