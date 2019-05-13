classdef Generator < handle
    % Class to to construct a time-dependent generator for a ctmc
    
    properties
        num_states
        num_steps
        time
        time_step
        rows            % contains the row indices of non-zero generator elements
        columns         % contains the column inidces of non-zero elements
        values          % each column contains the values of the generator corresponding to the time step 
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Generator(time,generator)
            % Store properties of the system and construct keymap and
            % generator for numerical integration
            obj.time = time;
            obj.num_states = size(generator,1);
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
            obj.values = repmat(value,[1,length(obj.time)]);
        end
        
        %% forward integration
        
        % integrate for given initial dist
        function [dist,t] = integrate(obj,initial_dist)
            ode_fun = @(t,y) obj.master_equation(t,y);
            [t,dist] = ode45(ode_fun,obj.time,initial_dist);
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
        function Q = interp_generator(obj,t)
            % get the next smallest time
            ind = floor(t/obj.time_step)+1;
            t_lower = obj.time_step*(ind-1);
            % interpolte the generator linearly
            alpha = (t-t_lower)/obj.time_step;
            Q = (1-alpha)*obj.values(:,ind)+alpha*obj.values(:,ind+1);
        end
        
    end
end

