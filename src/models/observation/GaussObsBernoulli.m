classdef GaussObsBernoulli < handle
    % Special Gaussian observation model for the intensity of the stemloops
    
    properties
        sigma
        num_nodes
        encoding
    end
    
    methods
        
        % constructor
        function obj = GaussObsBernoulli(sigma,num_nodes,encoding)
            obj.sigma = sigma;
            obj.num_nodes = num_nodes;
            obj.encoding = encoding;
        end
        
        % compute terminal condition of the backward equation
        function terminal = get_terminal(obj,moments,observation,sample_time)
            % distinguish between full and reduced state representation
            if strcmp(obj.encoding,'full')
                % compute sumary stats of the process
                exp_val = moments(2:2:end);
                % get update from first order moments
                terminal = (observation-exp_val)/obj.sigma^2;
                terminal = terminal-0.5/obj.sigma^2*(1-2*exp_val);
                % expand two binary
                terminal = [zeros(size(terminal)),terminal]';
                terminal = terminal(:);
            elseif strcmp(obj.encoding,'reduced')
                % get update from first order moments
                terminal = (observation-moments)/obj.sigma^2;
                terminal = terminal-0.5/obj.sigma^2*(1-2*moments);
            end
        end
        
        % compute the residual
        function residual = get_residual(obj,moments,observation,sample_time)
            % compute sumary stats of the process
            if strcmp(obj.encoding,'full')
                exp_val = moments(2:2:end);
            elseif strcmp(obj.encoding,'reduced')
                exp_val = moments;
            end
            % contribution of the normalizer
            residual = obj.num_nodes*0.5*log(2*pi)+obj.num_nodes*log(obj.sigma);
            % central moment contribution
            residual = residual+0.5/obj.sigma^2*sum(exp_val.*(1-exp_val));
            % observation contribution
            residual = residual+0.5/obj.sigma^2*sum((observation-exp_val).^2);
        end
        
        
        
        %         % evaluate log likilihood for a sample
        %         function llh = eval(obj,time,state,observation)
        %             % contribution of normalizer
        %             llh =  -0.5*log(2*pi)-log(obj.sigma);
        %             % contribution of observation
        %             llh = llh-0.5/obj.sigma^2*(observation-state(obj.species_selection))^2;
        %         end
        
    end
end

