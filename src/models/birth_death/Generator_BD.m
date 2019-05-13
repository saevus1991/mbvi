classdef Generator_BD < handle
    % Generator of birth death process with state-space truncation
    % disease model

    
    properties
        num_states
        system
        keymap
        generator
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Generator_BD(system)
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.system = system;
            obj.num_states = system.max_state+1;
            obj.keymap = obj.construct_keymap;
            obj.generator = obj.construct_generator;
        end
        
        % construct a keymap that associates each state of the graph with alabel
        function keymap = construct_keymap(obj)
            keymap = struct;
            keymap.index = [1:obj.num_states]';
            keymap.states = keymap.index-1;
        end
            
        % construct the generator in a sparse form   
        function generator = construct_generator(obj)
            % preparations
            row_index = [];
            col_index = [];
            values = [];
            %
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
            if state == 0
                prop = [obj.system.rates(1);0];
            elseif state == obj.system.max_state
                prop = [0;obj.system.rates(2)*obj.num_states];
            else
                prop = [obj.system.rates(1);obj.system.rates(2)*state];
            end
            % each nonzero object of prop corresponds to a possible transition
            target_ind = zeros(size(prop));
            for i = 1:length(prop)
                % get target states for nonzero props
                if prop(i) > 0
                    new_state = state+3-2*i;        % produces a +1 for birth reaction and a -1 for death reaction
                    target_ind(i) = obj.state2ind(new_state);
                end
            end
            % kill zeros from output
            target_prop = [prop(prop>0);-sum(prop)];
            target_ind = [target_ind(target_ind>0);ind];
            state_ind = ind*ones(size(target_prop));
        end
        
        %% simple helper functions
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.keymap.states(ind,:)';
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            ind = state+1;
        end
        
        % compute deterministic initial distribution
        function dist = initial_dist(obj,state)
            ind = obj.state2ind(state);
            dist = zeros(obj.num_states,1);
            dist(ind) = 1;
        end
        
        %% integrator
        
        % integrate an initial distribution forward in time using the generator
        function [t,dist] = integrate(obj,initial_dist,t_min,t_max)
            % set up ode function
            ode_fun = @(t,y) obj.generator*y;
            % integrate
            [t,dist] = ode45(ode_fun,[t_min,t_max],initial_dist);
            dist = dist';
        end
        
        
        %% summary statistics
        
        % compute the expected number of individuals
        function population_mean = population_mean(obj,dist)
            % initialize
            states = obj.ind2state(1:obj.num_states); 
            population_mean = sum(states'.*dist);
        end
        
    end
end

