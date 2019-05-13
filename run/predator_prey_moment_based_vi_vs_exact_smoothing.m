% load required libraries
addpath('../data/synthetic_data/')
addpath('../data/output')
addpath('../src/models/predator_prey/')

%% Description
% Produce plot from a predator prey model for the moment based vi

%% Load dataset

% load the simulations dataset
file_path = 'predator_prey_gauss_N_100.mat';
sim_data = load(file_path);
sim_data = sim_data.data_set;

% load the smoothed dataset
file_path = 'predator_prey_gauss_N_100_exact_smoothed.mat';
smoothed_data = load(file_path);
smoothed_data = smoothed_data.output;

% load the vi dataset
file_path = 'predator_prey_gauss_N_100_processed.mat';
mmvi_data = load(file_path);
mmvi_data = mmvi_data.output;

% number of files
N = length(sim_data.measurement);

%% load the model

system = sim_data.system;

% Gaussian observation noise
sigma = sim_data.sigma;
noise = sim_data.noise_model;

% save errors
error = zeros(N,2);

% pick file
for k = 1:100

% convert the moments
m = mmvi_data.moments{k};
m_t = mmvi_data.time{k}';
prey_mean_smooth =  m(:,1)';
prey_err_smooth = sqrt(m(:,3))';
pred_mean_smooth =  m(:,2)';
pred_err_smooth = sqrt(m(:,5))';
delta_t = m_t(2:end)-m_t(1:end-1);

% load corresponding moment equations of the smoother
smoother_mean = smoothed_data.smoother_stats{k,1};
smoother_time = smoothed_data.filter_time{k};
smoother_time(end) = smoother_time(end)+1e-3*(smoother_time(end)-smoother_time(end-1));
smoother_mean = discretize_data(smoother_time,smoother_mean,m_t);

% compute squared deviations
prey_error = sqrt((prey_mean_smooth-smoother_mean(1,:)).^2);
prey_error = 0.5*sum(delta_t.*(prey_error(2:end)+prey_error(1:end-1)))/(system.t_max-system.t_min);

pred_error = sqrt((pred_mean_smooth-smoother_mean(2,:)).^2);
pred_error = 0.5*sum(delta_t.*(pred_error(2:end)+pred_error(1:end-1)))/(system.t_max-system.t_min);

% store
error(k,:) = [prey_error,pred_error];

figure
hold on
plot(m_t,prey_mean_smooth)
plot(m_t,smoother_mean(1,:))
plot(m_t,pred_mean_smooth)
plot(m_t,smoother_mean(2,:))



end

