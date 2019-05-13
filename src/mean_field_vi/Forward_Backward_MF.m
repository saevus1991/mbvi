classdef Forward_Backward_MF < handle
    % Contains the generator of a CTMC in sparse form a well as keymaps and
    % the parameters structure
    
    properties
        num_states
        model
        rates
        rows                % contains the row indices of non-zero generator elements
        columns             % contains the column inidces of non-zero elements
        rate_map            % associates a rate constant to each row, column index pari
        blank_generator
        forward_generator
        backward_generator
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Forward_Backward_MF(model,rates)
            % Stor properties of the parameters and construct keymap and
            % generator for numerical integration
            obj.model = model;
            obj.rates = rates;
            obj.num_states = model.num_states;
            obj.construct_forward_generator();
        end
            
%         % construct the generator in a sparse form   
%         function construct_forward_generator(obj)
%             % preparations
%             row_index = [];
%             col_index = [];
%             values = [];
%             % iterate over all states
%             for i = 1:obj.num_states
%                 % get target states and propensities
%                 [target_ind,target_prop,state_ind,rate_ind] = obj.get_targets(i);
%                 % save in store
%                 row_index = [row_index;target_ind];
%                 col_index = [col_index;state_ind];
%                 values = [values;target_prop];
%             end
%             % set up generator
%             obj.rows = row_index;
%             obj.columns = col_index;
%             obj.forward_generator = values;
%         end
        
        % construct the generator as a value matrix from time_dependent rates   
        function construct_forward_generator(obj)
            % preparations
            row_index = [];
            col_index = [];
            rate_index = [];
            values = [];
            % iterate over all states
            for i = 1:obj.num_states
                % get target states and propensities
                [target_ind,target_prop,state_ind,rate_ind] = obj.get_targets(i);
                % save in store
                row_index = [row_index;target_ind];
                col_index = [col_index;state_ind];
                rate_index = [rate_index;rate_ind];
                values = [values;target_prop];
            end
            % set up generator
            obj.rows = row_index;
            obj.columns = col_index;
            obj.rate_map = rate_index;
            obj.blank_generator = values;
            obj.forward_generator = obj.update_generator(obj.rates);
        end
        
        % construct the generator in a sparse form   
        function values = update_generator(obj,rates)
            % preparations
            values = obj.blank_generator.*rates(obj.rate_map);
            for i = 1:obj.num_states
                values(obj.columns==i&obj.rows==i) = -sum(values(obj.columns==i&obj.rows~=i));
            end
        end
        
        % set up time varying generator    
        function values = timde_dep_generator(obj,mf_stats)
            % preparations
            values = obj.blank_generator.*obj.rates(obj.rate_map);
            values = repmat(values,[1,size(mf_stats{1},2)]);
            modified_values = values.*mf_stats{1}(obj.rate_map);
            values = values.*mf_stats{2}(obj.rate_map);
            for i = 1:obj.num_states
                values(obj.columns==i&obj.rows==i,:) = -sum(modified_values(obj.columns==i&obj.rows~=i),:);
            end
        end
        
        % get target state for state with index ind
        function [target_ind,target_prop,state_ind,rate_ind] = get_targets(obj,ind)
            % get state corresponding to index
            state = obj.ind2state(ind);
            % calculate propensities
            prop = obj.raw_propensity(state);
            % each nonzero object of prop corresponds to a possible transition
            target_ind = zeros(size(prop));
            for i = 1:length(prop)
                % get target states for nonzero props
                if prop(i) > 0
                    new_state = obj.update_state(state,i);
                    target_ind(i) = obj.state2ind(new_state);
                end
            end
            % index to associate rate constant
            rate_ind = (1:length(prop))';
            rate_ind = [rate_ind(prop>0);1];
            % kill zeros from output
            target_prop = [prop(prop>0);-sum(prop)];
            target_ind = [target_ind(target_ind>0);ind];
            state_ind = ind*ones(size(target_prop));
        end
        
        %% simple helper functions
        
        % compute propensity
        function prop = propensity(obj,state)
            prop = obj.model.propensity(state,obj.rates);
        end
        
        % compute raw propensity without rates
        function prop = raw_propensity(obj,state)
            prop = obj.model.raw_propensity(state);
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
        
        % integrate for given initial dist
        function dist = integrate(obj,initial_dist,tspan,values)
            ode_fun = @(t,y) obj.master_equation(t,y,tspan,values);
            [~,dist] = ode45(ode_fun,tspan,initial_dist);
            dist = dist';
        end
        
        % evaluates the r.h.s. of the master equation 
        function dydt = master_equation(obj,t,y,tspan,values)
            % get time-dependent generator
            Q = obj.interp_generator(t,tspan,values);
            % evaluate derivative
            dydt = sparse(obj.rows,obj.columns,Q,obj.num_states,obj.num_states)*y;
        end
        
        % interpolate the generator
        function Q = interp_generator(obj,t,tspan,values)
            % get the next smallest time
            time_step = tspan(2)-tspan(1);
            ind = floor((t-tspan(1))/time_step)+1;
            t_lower = time_step*(ind-1);
            % interpolte the generator linearly
            alpha = (t-t_lower)/time_step;
            try
            Q = (1-alpha)*values(:,ind)+alpha*values(:,ind+1);
            catch
                x = 1;
            end
        end
        
        % integrate for given initial dist
        function dist = backward_integrate(obj,terminal_dist,tspan,values)
            ode_fun = @(t,y) obj.backward_equation(t,y,tspan,values);
            [~,dist] = ode45(ode_fun,tspan(end:-1:1),terminal_dist);
            dist = dist';
            dist = dist(:,end:-1:1);
        end
        
        % evaluates the r.h.s. of the master equation 
        function dydt = backward_equation(obj,t,y,tspan,values)
            % get time-dependent generator
            Q = obj.interp_generator(t,tspan,values);
            % evaluate derivative
            dydt = -sparse(obj.rows,obj.columns,Q,obj.num_states,obj.num_states)'*y;
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

