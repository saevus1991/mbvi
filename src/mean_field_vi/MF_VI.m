classdef variational_engine < handle
    % Moment based variational inference for MJPs
    %   This version uses a fixed time grid
    
    properties
        model               % object that contains equations for forward, backward and gradient
        obs_model           % object that describes the observation model
        options 
        rates               % the sufficient statistics for the parameters 
        var_rates           % the hyperparameters of the variational distribution
        initial_time        % initial time for calculations
        final_time          % final time of observations
        time_step
        sample_times        % times at which the system is observed
        observations        % observed values at the sample time
        num_intervals       % number of intervals that have to be calculated
        time_grid           % the time grid on which all functions are evaluated and stored
        grid_size           % number of time points for each subinteral
        moments             % cell array for the momens and corresponding time grid
        propensities
        constraints         % cell array for the constraints and correpsonding time grid
        statistics          % cell containing the sufficient statistics for updating the parameters
        residual            % the <log(p(y|x)>
    end
    
    methods
        
        %% setup functions
        
        % constructor method
        function obj = variational_engine(model,obs_model,time_step,sample_times,observations,options)
            % pass on the values to member variables
            obj.model = model;
            obj.obs_model = obs_model;
            obj.options = options;
            obj.initial_time = model.tspan(1);
            obj.final_time = model.tspan(2);
            obj.time_step = time_step;
            obj.sample_times = sample_times;
            obj.observations = observations;
            % set number of time intervals
            obj.num_intervals = length(sample_times)+1;
            % construct the cell arrays
            obj.time_grid = cell(obj.num_intervals,1);
            obj.moments = cell(obj.num_intervals,obj.model.num_species);
            obj.constraints = cell(obj.num_intervals,obj.model.num_species);
            obj.statistics = cell(3,1);
            obj.residual = zeros(size(observations));
            % initialize the rates
            obj.initialize_rates();
            % construct time grid
            obj.set_time_grid()
            % initialzie the residuals 
            obj.evaluate_residuals();
            % initialize the statistics
            obj.initialize_statistics();
        end
        
        % initialize the rates from the prior
        function initialize_rates(obj)
            obj.rates = zeros(length(obj.model.initial_rates),2);
            % if mode is set to smoothing, initialize with fixed rates
            if strcmp(obj.options.operation_mode,'smoothing')
                obj.rates(:,1) = obj.model.initial_rates;
                obj.rates(:,2) = obj.model.initial_rates;
                % if mode is set to inference, use initilization that accounts for expectation
            elseif strcmp(obj.options.operation_mode,'inference')
                obj.rates(:,1) = obj.model.alpha./obj.model.beta;
                obj.rates(:,2) = exp(psi(obj.model.alpha))./obj.model.beta;
                % also initialize rates of the variational hyper parameters
                obj.var_rates = zeros(length(obj.model.alpha),2);
                obj.var_rates(:,1) = obj.model.alpha;
                obj.var_rates(:,2) = obj.model.beta;
            end
        end
        
        % set up the time grid
        function set_time_grid(obj,time_step)
            if nargin == 1
                time_step = obj.time_step;
            end
            % initialize the first interval
            num_steps = ceil((obj.sample_times(1)-obj.initial_time)/time_step);
            obj.time_grid{1} = linspace(obj.initial_time,obj.sample_times(1),num_steps)';
            % iterate over intermediate intervals
            for i = 2:obj.num_intervals-1
                num_steps = ceil((obj.sample_times(i)-obj.sample_times(i-1))/time_step);
                obj.time_grid{i} = linspace(obj.sample_times(i-1),obj.sample_times(i),num_steps)';
            end
            % initialize last interval
            num_steps = ceil((obj.final_time-obj.sample_times(end))/time_step);
            obj.time_grid{obj.num_intervals} = linspace(obj.sample_times(obj.num_intervals-1),obj.final_time,num_steps)';
            % compute number of time points per interval
            obj.grid_size = cell2mat(cellfun(@length,obj.time_grid,'un',0));
            obj.modify_time_grid();
        end
        
        % modify time grid
        function modify_time_grid(obj)
            % slightly disturb the time step at the end of the grid values
            % such that it is different from the first 
            for i = 1:obj.num_intervals
                obj.time_step = obj.time_grid{i}(2)-obj.time_grid{i}(1);
                obj.time_grid{i}(end) = obj.time_grid{i}(end)-1e-6*obj.time_step;
            end
        end
        
        % initialize statistics
        function initialize_statistics(obj)
            obj.statistics{1} = zeros(obj.num_intervals,obj.model.num_reactions);
            obj.statistics{2} = zeros(obj.num_intervals,obj.model.num_reactions);
            obj.statistics{3} = zeros(obj.num_intervals,obj.model.num_reactions);
            obj.compute_statistics();
        end
        
        %% basic update functions
        
        % forward update
        function species_update(obj,ind)
            % preparations
            initial = obj.model.initial{ind};
            % iterate over time intervals for forward solution
            for i = 1:obj.num_intervals
                forward_update = @(t,y) obj.model.moments(t,y,obj.time_grid{i},obj.constraints{i,ind});
                [~,moment] = ode45(forward_update,obj.time_grid{i}',initial);
                obj.moments{i,ind} = moment;
                initial = obj.moments{i,ind}(end,:)';
            end
            % set constraint in last interval to zero
            num_constraints = length(obj.model.initial{ind});
            obj.constraints{obj.num_intervals,ind} = zeros(length(obj.time_grid{obj.num_intervals}),num_constraints);
            final = zeros(size(obj.model.initial));
            % iterate backward
            for i = obj.num_intervals-1:-1:1
                backward_update = @(t,y)obj.model.constraints(t,y,obj.time_grid{i},obj.moments{i},obj.rates);
                final = final+obj.obs_model.get_terminal(obj.moments{i}(end,:)',obj.observations(:,i),obj.sample_times(i));
                [~,constraint] = ode45(backward_update,obj.time_grid{i}(end:-1:1)',final);
                obj.constraints{i,ind} = constraint(end:-1:1,:);
                final = obj.constraints{i,ind}(1,:)';
            end            
        end
        
        % compute statistics
        function compute_statistics(obj)
            % iterate over the intervals to compute the sufficient statistics
            for i = 1:obj.num_intervals
                % evaluate sumary statistics
                delta_t = obj.time_grid{i}(2:end)-obj.time_grid{i}(1:end-1);
                stat1 = obj.rates(:,1)'.*obj.propensities{i};
                stat1 = 0.5*(stat1(2:end,:)+stat1(1:end-1,:));
                stat2 = obj.control{i}.*obj.propensities{i};
                stat2 = 0.5*(stat2(2:end,:)+stat2(1:end-1,:));
                stat3 = obj.control{i}.*obj.propensities{i}.*log(obj.control{i}./obj.rates(:,2)');
                stat3 = 0.5*(stat3(2:end,:)+stat3(1:end-1,:));
                % remove nans
                stat3(isnan(stat3)) = 0;
                % store
                obj.statistics{1}(i,:) = sum(stat1.*delta_t);
                obj.statistics{2}(i,:) = sum(stat2.*delta_t);
                obj.statistics{3}(i,:) = sum(stat3.*delta_t);
            end                        
        end
        
        % get the sufficient statistics
        function [stat1,stat2] = get_statistics(obj)
            stat1 = sum(obj.statistics{1});
            stat2 = sum(obj.statistics{2});
        end
        
        % evalute the residuals for the gaussian noise model and currently
        % only for one dimensional problems
        function evaluate_residuals(obj)
            res = zeros(1,size(obj.observations,2));
            for i = 1:obj.num_intervals-1
                res(:,i) = obj.obs_model.get_residual(obj.moments{i}(end,:)',obj.observations(:,i),obj.sample_times(i));
            end
            obj.residual = res;
        end
        
        % evaluate the evidence lower bound
        function [kl] = evaluate_kl(obj)            
            % compute sufficient statistics 
            stat1 = sum(obj.statistics{1});
            stat2 = sum(obj.statistics{2});
            stat3 = sum(obj.statistics{3});
            % integral contribution
            kl = sum(stat1-stat2+stat3);
            % residual contribution
            kl = kl+sum(obj.residual(:));
        end
        
               
        %% objective functions
        
        % target function for optimization
        function [elbo,grad] = objective_function(obj,control)
            % set control
            if nargin == 2
                obj.control = mat2cell(control,obj.grid_size);
            end
            % evaluate backward integration only if gradient information is required 
            if nargout == 2
                % update all sub-components
                obj.moment_update();
                obj.constraint_update();
                obj.gradient_update();
                % compute the elbo
                obj.evaluate_residuals();
                obj.compute_statistics();
                elbo = obj.evaluate_kl;
                % transform the gradient to the required form
                grad = cell2mat(obj.gradient);
            elseif nargout == 1
                if nargin == 2 && any(control(:)<0)
                    elbo = Inf;
                else
                    % update only the required sub-components
                    obj.moment_update();
                    % compute the elbo
                    obj.evaluate_residuals();
                    obj.compute_statistics();
                    elbo = obj.evaluate_kl;
                end
            end
        end
        
        function [elbo] = full_objective_function(obj,var_par,contr) 
            % set new control
            obj.control = mat2cell(contr,obj.grid_size);
            % set new variational parameters and the rates
            obj.var_rates = var_par;
            obj.rates(:,1) = var_par(:,1)./var_par(:,2);
            obj.rates(:,2) = exp(psi(var_par(:,1)))./var_par(:,2);
            % update the dynamic equations
            obj.moment_update();
            obj.constraint_update();
            obj.gradient_update();
            % compute kl objective function
            obj.evaluate_residuals();
            obj.compute_statistics();
            elbo = obj.kl_gamma(var_par)+obj.evaluate_kl;
        end
        
        %% simple getters    
        
        % get the moments
        function [m_t,m] = get_moments(obj)
            m_t = [];
            m = [];
            for i = 1:obj.num_intervals
                m_t = [m_t;obj.time_grid{i}];
                m = [m;obj.moments{i}];
            end
        end
        
        % get the constraints
        function [c_t,c] = get_constraints(obj)
            c_t = [];
            c = [];
            for i = 1:obj.num_intervals
                c_t = [c_t;obj.time_grid{i}];
                c = [c;obj.constraints{i}];
            end
        end
        
        % get the control;
        function [a_t,a] = get_control(obj)
            a_t = [];
            a = [];
            for i = 1:obj.num_intervals
                a_t = [a_t;obj.time_grid{i}];
                a = [a;obj.control{i}];
            end
        end
        
        % get the gradient
        function [g_t,g] = get_gradient(obj)
            g_t = [];
            g = [];
            for i = 1:obj.num_intervals
                g_t = [g_t;obj.time_grid{i}];
                g = [g;obj.gradient{i}];
            end
        end
        
        % update the rates using the sufficient statistics
        function [stat1,stat2] = get_stats(obj)
            % compute sufficient statistics 
            stat1 = sum(obj.statistics{1});
            stat2 = sum(obj.statistics{2});
        end
        
    end
end

