classdef Generator_GE < handle
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
        function  obj = Generator_GE(max_state)
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.max_state = max_state;
            obj.num_states = prod(max_state+1);
            obj.keymap = obj.construct_keymap;
        end
        
        % construct a keymap that associates each state with an index
        function keymap = construct_keymap(obj)
            keymap = struct;
            keymap.index = [1:obj.num_states]';
            states = combvec(0:obj.max_state(3),0:obj.max_state(2),0:obj.max_state(1))';
            keymap.states = states(:,end:-1:1);
        end          
        
        %% simple helper functions
        
        % evalute propensity
        function prop = propensity(obj,state,rates)
            % initialize
            prop = zeros(size(rates));
            % normal mass action 
            prop(1) = (1-state(1));
            prop(2) = state(1);
            prop(3) = state(1);
            prop(4) = state(2);
            prop(5) = state(2);
            prop(6) = state(3);
            prop = prop.*rates;
            % correct if state hits max values
            if state(2) == obj.max_state(2)
                prop(3) = 0;
            end
            if state(3) == obj.max_state(3)
                prop(5) = 0;
            end
        end
        
        % compute target state
        function new_state = update_state(obj,state,index)
            new_state = state;
            switch index
                case 1
                    new_state(1) = 1;
                case 2
                    new_state(1) = 0;
                case 3
                    new_state(2) = new_state(2)+1;
                case 4
                    new_state(2) = new_state(2)-1;
                case 5
                    new_state(3) = new_state(3)+1;
                case 6
                    new_state(3) = new_state(3)-1;
            end
        end
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.keymap.states(ind,:)';
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            %[~,ind] = ismember(state',obj.keymap.states,'rows');
            ind = state(1)*(obj.max_state(2)+1)*(obj.max_state(3)+1)+state(2)*(obj.max_state(3)+1)+state(3)+1;
        end
        
        
        %% summary statistics
        
        % compute the expected number of individuals
        function population_mean = population_mean(obj,dist)
            % initialize
            states = obj.ind2state(1:obj.num_states)'; 
            gene_mean = sum(states(:,1).*dist);
            mrna_mean = sum(states(:,2).*dist);
            protein_mean = sum(states(:,3).*dist);
            population_mean = [gene_mean;mrna_mean,protein_mean];
        end        
        
    end
end


       
