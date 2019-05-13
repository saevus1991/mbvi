classdef Forward_Backward < handle
    % Contains the generator of a CTMC in sparse form a well as keymaps and
    % the parameters structure
    
    properties
        num_states
        model
        parameters
        generator
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Forward_Backward(model,parameters)
            % Stor properties of the parameters and construct keymap and
            % generator for numerical integration
            obj.model = model;
            obj.parameters = parameters;
            obj.num_states = model.num_states;
            obj.generator = obj.construct_generator;
        end
            
        % construct the generator in a sparse form   
        function generator = construct_generator(obj,rates)
            if nargin == 2
                obj.parameters.rates = rates;
            end
            % preparations
            row_index = [];
            col_index = [];
            values = [];
            % iterate over all states
            for i = 1:obj.num_states
                % get target states and propensities
                [target_ind,target_prop,state_ind] = obj.get_targets(i);
                % save in store
                row_index = [row_index;target_ind];
                col_index = [col_index;state_ind];
                values = [values;target_prop];
            end
            % set up generator
            generator = sparse(row_index,col_index,values,obj.num_states,obj.num_states);
        end
        
        % get target state for state with index ind
        function [target_ind,target_prop,state_ind] = get_targets(obj,ind)
            % get state corresponding to index
            state = obj.ind2state(ind);
            % calculate propensities
            prop = obj.propensity(state);
            % each nonzero object of prop corresponds to a possible transition
            target_ind = zeros(size(prop));
            for i = 1:length(prop)
                % get target states for nonzero props
                if prop(i) > 0
                    new_state = obj.update_state(state,i);
                    target_ind(i) = obj.state2ind(new_state);
                end
            end
            % kill zeros from output
            target_prop = [prop(prop>0);-sum(prop)];
            target_ind = [target_ind(target_ind>0);ind];
            state_ind = ind*ones(size(target_prop));
        end
        
        %% simple helper functions
        
        % compute propensity
        function prop = propensity(obj,state)
            prop = obj.model.propensity(state,obj.parameters.rates);
        end
        
        % update state if reaction revent event happens
        function new_state = update_state(obj,state,event)
            new_state = obj.model.update_state(state,event);
        end
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.model.ind2state(ind);
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            ind = obj.model.state2ind(state);
        end
        
        % compute deterministic initial distribution
        function dist = initial_dist(obj,state)
            ind = obj.state2ind(state);
            dist = zeros(obj.num_states,1);
            dist(ind) = 1;
        end
        
        %% integrator
        
        % integrate an initial distribution forward in time using the generator
        function [t,dist] = integrate(obj,initial_dist,t_span)
            % set up ode function
            ode_fun = @(t,y) obj.generator*y;
            % integrate
            [t,dist] = ode45(ode_fun,t_span,initial_dist);
            dist = dist';
        end
        
        % integrate a distribution backward in time using the backward generator
        function [t,dist] = backward_integrate(obj,terminal_dist,t_span)
            % set up ode function
            ode_fun = @(t,y) -obj.generator'*y;
            % integrate
            [t,dist] = ode45(ode_fun,t_span,terminal_dist);
            dist = dist';
        end        
        
        %% observation related
        
         % compute observation likelihood for a given observation model
        function obs_likelihood = observation_likelihood(obj,obs_model,observation)
            % initialize
            states = obj.ind2state(1:obj.num_states); 
            obs_likelihood = zeros(obj.num_states,size(observation,2));
            % iteratore over observations
            for i = 1:size(observation,2)
                obs_likelihood(:,i) = obs_model(states,observation(:,i));
            end
        end       

    end
end

