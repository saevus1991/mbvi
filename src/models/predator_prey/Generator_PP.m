classdef Generator_PP < handle
    % Generator of autocatalytic birth death process with state-space
    % truncation
    
    properties
        max_state
        num_states
        keymap
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor
        function  obj = Generator_PP(max_state)
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.max_state = max_state;
            obj.num_states = (max_state(1)+1)*(max_state(2)+1);
            obj.keymap = obj.construct_keymap;
        end
        
        % construct a keymap that associates each state with an index
        function keymap = construct_keymap(obj)
            keymap = struct;
            keymap.index = [1:obj.num_states]';
            prey = repelem((0:obj.max_state(1))',obj.max_state(2)+1);
            predator = repmat((0:obj.max_state(2))',[obj.max_state(1)+1,1]);
            keymap.states = [prey,predator];
        end          
        
        %% simple helper functions
        
        % evalute propensity
        function prop = propensity(obj,state,rates)
            % initialize
            prop = zeros(size(rates));
            % normal mass action 
            prop(1) = state(1);
            prop(2) = state(1)*state(2);
            prop(3) = state(1)*state(2);
            prop(4) = state(2);
            prop = prop.*rates;
            % correct if state hits max values
            if state(1) == obj.max_state(1)
                prop(1) = 0;
            end
            if state(2) == obj.max_state(2)
                prop(3) = 0;
            end
        end

        % evalute raw propensity without rates
        function prop = raw_propensity(obj,state)
            % initialize
            prop = zeros(4);
            % normal mass action 
            prop(1) = state(1);
            prop(2) = state(1)*state(2);
            prop(3) = state(1)*state(2);
            prop(4) = state(2);
            % correct if state hits max values
            if state(1) == obj.max_state(1)
                prop(1) = 0;
            end
            if state(2) == obj.max_state(2)
                prop(3) = 0;
            end
        end        

        % compute target state
        function new_state = update_state(obj,state,index)
            new_state = state;
            switch index
                case 1
                    new_state(1) = new_state(1)+1;
                case 2
                    new_state(1) = new_state(1)-1;
                case 3
                    new_state(2) = new_state(2)+1;
                case 4
                    new_state(2) = new_state(2)-1;
            end
        end
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.keymap.states(ind,:)';
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            ind = state(1,:)*(obj.max_state(2)+1)+1+state(2,:);
        end
        
        
        %% summary statistics
        
        % compute the expected number of individuals
        function population_mean = population_mean(obj,dist)
            % initialize
            states = obj.ind2state(1:obj.num_states)'; 
            prey_mean = sum(states(:,1).*dist);
            predator_mean = sum(states(:,2).*dist);
            population_mean = [prey_mean;predator_mean];
        end        
        
        % compute the expected number of individuals
        function population_std = population_std(obj,dist,population_mean)
            % initialize
            states = obj.ind2state(1:obj.num_states)'; 
            prey_std = sum((states(:,1)-population_mean(1,:)).^2.*dist);
            predator_std = sum((states(:,2)-population_mean(2,:)).^2.*dist);
            population_std = [prey_std;predator_std];
        end      
        
    end
end


       
