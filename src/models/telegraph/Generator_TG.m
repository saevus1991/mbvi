classdef Generator_TG < handle
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
        function  obj = Generator_TG()
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.max_state = 1;
            obj.num_states = max_state+1;
            obj.keymap = obj.construct_keymap;
        end
        
        % construct a keymap that associates each state with an index
        function keymap = construct_keymap(obj)
            keymap = struct;
            keymap.index = [1:obj.num_states]';
            keymap.states = keymap.index-1;
        end
               
        %% simple helper functions
        
        % compute the propensity
        function prop = propensity(obj,state,rates)
            if state == 0
                prop = rates(1);
            elseif state == 1
                prop = rates(2);
            end
        end

        % update the state
        function new_state = update_state(obj,state,event)
            switch event
                case 1
                    new_state = 1;
                case 2
                    new_state = 0;
            end
        end        
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.keymap.states(ind,:)';
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            ind = state+1;
        end
        
    end
end

