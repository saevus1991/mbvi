%% test ssa for SIS disease model on a graph

addpath('../src/models/autocatalytic_birth_death')
addpath('../src/exact_integration')
addpath('../src/ssa')

% specify graph by adjecency matrix

 % fix state space trunction
 max_state = 120;
 
 % construct model
 abd_model = Generator_ABD(max_state);
 
 % set up a parameter structure for the model
 system.initial = 5;
 system.rates = [0.45;0.4];      % infection rate, recovery rate
 system.t_min = 0;
 system.t_max = 200;
 
 % set up integrator object
 generator = Forward_Backward(abd_model,system);
 
 % get initial distribution
 initial_dist = generator.initial_dist(system.initial);
 
 % solve master equation
 [t,dist] = generator.integrate(initial_dist,[system.t_min,system.t_max]);
 
 % compute marginal
 population_mean = generator.population_mean(dist);
 
 % set up synthetic observations
 obs_times = [20,30,40,50,60,80];
 obs_vals = [25,40,60,62,55,62];
 
 % gaussian observation model
 sigma = 5;
 obs_model = @(x,y) -0.5*log(2*pi*sigma^2)-0.5*(y-x).^2/sigma^2;
 obs_likelihood = generator.observation_likelihood(obs_model,obs_vals);
 filter_norm = zeros(1,size(obs_likelihood,2));
 
 % construct sampling times
 delta_t = 0.1;
 sample_times = cell(length(obs_times)+1,1);
 sample_times{1} = system.t_min:delta_t:(obs_times(1)-0.001*delta_t);
 for i = 2:length(sample_times)-1
     sample_times{i} = obs_times(i-1):delta_t:(obs_times(i)-0.001*delta_t);
 end
 sample_times{end} = obs_times(end):delta_t:(system.t_max-0.001*delta_t);
 
 % compute filtering distribution
 forward_dist = cell(size(sample_times));
 for i = 1:length(sample_times)
     [~,dist] = generator.integrate(initial_dist,sample_times{i});
     forward_dist{i} = dist;
     if i < length(sample_times)
        initial_dist = obs_likelihood(:,i)+log(dist(:,end));
        max_log = max(initial_dist);
        initial_dist = exp(initial_dist-max_log);
        norm = sum(initial_dist);
        initial_dist = initial_dist/norm;
        filter_norm(i) = max_log+log(norm);
     end
 end
 
 % compute backward filter
 backward_dist = cell(size(sample_times));
 backward_dist{end} = ones(length(initial_dist),size(sample_times{end},2));
 terminal_dist = exp(obs_likelihood(:,i-1)-filter_norm(i-1));
 for i = (length(sample_times)-1):-1:1
     [~,dist] = generator.backward_integrate(terminal_dist,sample_times{i}(end:-1:1));
     backward_dist{i} = dist(:,end:-1:1);
     if i > 1
        terminal_dist = log(dist(:,end))+obs_likelihood(:,i-1)-filter_norm(i-1);
        terminal_dist = exp(terminal_dist);
     end
 end
 
 % compute the mean of the filtering distribution
 filter_time = cell2mat(sample_times');
 filter_mean = generator.population_mean(cell2mat(forward_dist'));
 
 % compute the mean of the smoothing disribution
 smoother_mean = generator.population_mean(cell2mat(forward_dist').*cell2mat(backward_dist'));
 
 % simulate trajectory from sparse gillespie implementation
 [sim_times,sim_states] = ssa_sparse(system.initial,[system.t_min,system.t_max],generator.generator);
 sim_states = generator.ind2state(sim_states);
 sim_states = discretize_data([system.t_min,sim_times],[system.initial,sim_states],t);
 
 % plot marginals
 figure
 hold on
 plot(t,population_mean,'-k')
 plot(filter_time,filter_mean,'-b')
 plot(filter_time,smoother_mean,'-r')
 plot(t,sim_states,'-g')
 plot(obs_times,obs_vals,'kx')
 xlabel('Time','fontsize',14)
 

 
 % simulate example trajectory
 
 %[times,states] = ssa_sis(system);
 
