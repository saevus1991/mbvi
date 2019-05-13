classdef Constraint < handle
    % Class to to construct a time-dependent constraints that have to be solved backward in time
    
    properties
        num_states
        num_steps
        time
        time_step
        rows            % contains the row indices of non-zero generator elements
        columns         % contains the column inidces of non-zero elements
        values_exp          % each column contains the values of the generator corresponding to the time step 
        values_log_exp
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Constraint(time,generator)
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.time = time;
            obj.num_states = num_states;
            obj.num_steps = length(time);
            obj.time_step = time(2)-time(1);
            % initialize with static generator
            obj.initialize_generator(generator);
        end
        
        % initialize the time dependent generator from a constant initial generator provided in sparse martrix form
        function initialize_generator(obj,generator)
            % extract row form
            [row,column,value] = find(generator);
            % store in corresponding object properties
            obj.rows = row;
            obj.columns = column;
            obj.values = zeros(length(value),length(obj.time));
        end
        
        %% forward integration
        
        % integrate for given initial dist
        function [dist,t] = integrate(obj,initial_dist)
            [t,dist] = ode45(obj.master_equation,obj.time,initial_dist);
            dist = dist';
        end
        
        % evaluates the r.h.s. of the master equation 
        function dydt = master_equation(obj,t,y)
            % get time-dependent generator
            Q = obj.interp_generator(t);
            % evaluate derivative
            dydt = sparse(obj.rows,obj.columns,Q,obj.num_states,obj.num_states)*y;
        end
        
        % interpolate the generator
        function Q = interp_generator(t)
            % get the next smallest time
            ind = floor(t/delta_t)+1;
            t_lower = delta_t*(ind-1);
            % interpolte the generator linearly
            alpha = (t-t_lower)/delta_t;
            Q = (1-alpha)*obj.values(:,ind)+alpha*obj.values(:,ind+1);
        end
        
    end
end

