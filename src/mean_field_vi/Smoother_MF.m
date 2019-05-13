classdef Smoother_MF < handle
    % Contains the generator of a CTMC in sparse form a well as keymaps and
    % the parameters structure
    
    properties
        initial_time        % initial time for calculations
        final_time          % final time of observations
        sample_times        % times at which the system is observed
        observations        % observed values at the sample time
        filter_norm         % stores normalization constant of forward filtering
        num_intervals       % number of intervals that have to be calculated
        time_grid           % the time grid on which all functions are evaluated and stored
        grid_size           % number of time points for each subinteral
        forward_dist        % forward filtered distribution 
        backward_dist       % result of backward filtering
        forward_generator   % time dependent forward generator
        backward_generator  % time depenent backward generator
        observation_likelihood      % precomputed observation likelihood required to reset filters at observation times
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor method
        function obj = Smoother_MF(initial_time,final_time,time_step,sample_times,observations)
            % pass on the values to member variables
            obj.initial_time = initial_time;
            obj.final_time = final_time;
            obj.sample_times = sample_times;
            obj.observations = observations;
            % set number of time intervals
            obj.num_intervals = length(sample_times)+1;
            % construct the cell arrays
            obj.time_grid = cell(obj.num_intervals,1);
            obj.forward_dist = cell(obj.num_intervals,1);
            obj.backward_dist = cell(obj.num_intervals,1);
            obj.forward_generator = cell(obj.num_intervals,1);
            obj.backward_generator = cell(obj.num_intervals,1);
            % construct time grid
            obj.set_time_grid(time_step)
        end
            
        % set up the time grid
        function set_time_grid(obj,time_step)
            % initialize the first interval
            num_steps = ceil((obj.sample_times(1)-obj.initial_time)/time_step);
            obj.time_grid{1} = linspace(obj.initial_time,obj.sample_times(1),num_steps);
            % iterate over intermediate intervals
            for i = 2:obj.num_intervals-1
                num_steps = ceil((obj.sample_times(i)-obj.sample_times(i-1))/time_step);
                obj.time_grid{i} = linspace(obj.sample_times(i-1),obj.sample_times(i),num_steps);
            end
            % initialize last interval
            num_steps = ceil((obj.final_time-obj.sample_times(end))/time_step);
            obj.time_grid{obj.num_intervals} = linspace(obj.sample_times(obj.num_intervals-1),obj.final_time,num_steps);
            % compute number of time points per interval
            obj.grid_size = cell2mat(cellfun(@length,obj.time_grid,'un',0));
            obj.modify_time_grid();
        end
        
        % modify time grid
        function modify_time_grid(obj)
            % slightly disturb the time step at the end of the grid values
            % such that it is different from the first 
            for i = 1:obj.num_intervals
                time_step = obj.time_grid{i}(2)-obj.time_grid{i}(1);
                obj.time_grid{i}(end) = obj.time_grid{i}(end)-1e-6*time_step;
            end
        end
        
        % initialize forward generator
        function initialize_forward_generator(obj,forward_backward)
            for i = 1:obj.num_intervals
                obj.forward_generator{i} = repmat(forward_backward.forward_generator,[1,size(obj.time_grid{i},2)]);
            end
        end
        
        %% getters 

        % get the time grid into single matrix
        function time = get_time(obj)
            % get time
            time = cell2mat(obj.time_grid');
        end        
        
        % get the filter distribution into single matrix
        function forward_dist = get_forward(obj)
            % get distribution
            forward_dist = cell2mat(obj.forward_dist');
        end
        
         % get the backward distribution into single matrix
        function backward_dist = get_backward(obj)
            % get distribution
            backward_dist = cell2mat(obj.backward_dist');
        end
        
         % get the backward distribution into single matrix
        function smoothed_dist = get_smoothed(obj,index)
            % get distribution
            if nargin == 1
                smoothed_dist = cell2mat(obj.backward_dist').*cell2mat(obj.forward_dist');
            elseif nargin == 2
                smoothed_dist = obj.backward_dist{index}.*obj.forward_dist{index};
            end
        end
        
        % get tehe log evidence estimate
        function evidence = get_evidence(obj)
            evidence = sum(obj.filter_norm);
        end
        
        %% setter
        
        % set observation likelihood
        function set_obs_likelihood(obj,model,obs_model)
            % initialize
            states = model.ind2state(1:model.num_states)'; 
            obs_likelihood = zeros(model.num_states,size(obj.observations,2));
            % iterate over observations
            for i = 1:size(obj.observations,2)
                obs_likelihood(:,i) = obs_model(states,obj.observations(:,i));
            end
            obj.observation_likelihood = obs_likelihood;
        end
        
        %% integrators
        
        % perform the forward filter
        function forward_filter(obj,forward_backward,initial_dist)
            % iterate over intervals defined by the time points of
            % observations
            for i = 1:obj.num_intervals
                dist = forward_backward.integrate(initial_dist,obj.time_grid{i},obj.forward_generator{i});
                obj.forward_dist{i} = dist;
                if i < obj.num_intervals
                    initial_dist = obj.observation_likelihood(:,i)+log(abs(dist(:,end)));
                    max_log = max(initial_dist);
                    initial_dist = exp(initial_dist-max_log);
                    norm = sum(initial_dist);
                    initial_dist = initial_dist/norm;
                    obj.filter_norm(i) = max_log+log(norm);
                end
            end
        end
        
        % perform the backward filter
        function backward_filter(obj,forward_backward)
            % the interval after the last observation becomes constant one
            obj.backward_dist{end} = ones(size(obj.observation_likelihood,1),size(obj.time_grid{end},2));
            terminal_dist = exp(obj.observation_likelihood(:,end)-obj.filter_norm(end));
            % iterate backward over the remaining intervals
            for i = (obj.num_intervals-1):-1:1
                dist = forward_backward.backward_integrate(terminal_dist,obj.time_grid{i},obj.backward_generator{i});
                obj.backward_dist{i} = dist;
                if i > 1
                    terminal_dist = log(abs(dist(:,1)))+obj.observation_likelihood(:,i-1)-obj.filter_norm(i-1);
                    terminal_dist = exp(terminal_dist);
                end
            end
        end
                
        
    end
end

