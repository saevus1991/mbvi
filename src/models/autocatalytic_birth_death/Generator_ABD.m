classdef Generator_ABD < handle
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
        function  obj = Generator_ABD(max_state)
            % Stor properties of the system and construct keymap and
            % generator for numerical integration
            obj.max_state = max_state;
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
                prop = [0;0];
            elseif state == obj.max_state
                prop = [0;rates(2)*obj.num_states];
            else
                prop = [rates(1);rates(2)]*state;
            end
        end
        
        % compute the propensity independent of thre rates and return an index of the corresponding rates
        function prop = raw_propensity(obj,state)
            if state == 0
                prop = [0;0];
            elseif state == obj.max_state
                prop = [0;obj.num_states];
            else
                prop = [1;1]*state;
            end
        end

        % update the state
        function new_state = update_state(obj,state,event)
            new_state = state+3-2*event;
        end        
        
        % convert index to state
        function state = ind2state(obj,ind)
            state = obj.keymap.states(ind,:)';
        end
        
        % convert state to index
        function ind = state2ind(obj,state)
            ind = state+1;
        end
        
        %% summary statistics
        
        % compute the expected number of individuals
        function population_mean = population_mean(obj,dist)
            % initialize
            states = obj.ind2state(1:obj.num_states)'; 
            population_mean = sum(states.*dist);
        end   
        
        % compute the expected number of individuals
        function [mf_stats1,mf_stats2] = mf_vi_stats(obj,dist,tol)
            % initialize
            states = obj.ind2state(1:obj.num_states)'; 
            states(1) = tol;
            mf_stats1 = sum(states.*dist);
            mf_stats2 = exp(sum(log(states).*dist));
        end   
        
    end
end

