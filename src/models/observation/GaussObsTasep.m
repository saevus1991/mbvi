classdef GaussObsTasep < handle
    % Special Gaussian observation model for the intensity of the stemloops
    
    properties
        obs_param
        alpha
    end
    
    methods
        
        % constructor
        function obj = GaussObsTasep(obs_param,alpha)
            obj.obs_param = obs_param;
            obj.alpha = alpha;
        end
        
        % compute terminal condition of the backward equation
        function terminal = get_terminal(obj,moments,observation,sample_time)
            % compute required sumary stats
            [I,grad_I] = obj.expected_intensity(moments,sample_time);
            [~,grad_Isq] = obj.variance_intensity(moments,sample_time);
            % get update from first order moments
            terminal = (observation-I)/obj.obs_param(5)^2*grad_I;
            terminal = terminal-0.5/obj.obs_param(5)^2*grad_Isq;
        end
        
        % compute the residual
        function residual = get_residual(obj,moments,observation,sample_time)
            % compute sumary stats of the process
            exp_intensity = obj.expected_intensity(moments,sample_time);
            var_intensity = obj.variance_intensity(moments,sample_time);
            % contribution of the normalizer
            residual = 0.5*log(2*pi)+log(obj.obs_param(5));
            % central moment contribution
            residual = residual+0.5/obj.obs_param(5)^2*var_intensity;
            % observation contribution
            residual = residual+0.5/obj.obs_param(5)^2*(observation-exp_intensity).^2;
        end
        
        %% helper function for calculatio
        
        function [I,grad_I] = expected_intensity(obj,moments,sample_time)
            % compute intensity from moments
            I = obj.obs_param(1)+exp(-obj.obs_param(3)*sample_time).*(obj.obs_param(2)+obj.obs_param(4)*obj.alpha'*moments);
            % if required, evaluate gradient
            if nargout == 2
                grad_I = exp(-obj.obs_param(3)*sample_time)*obj.obs_param(4)*obj.alpha;
            end
        end
        
        function [Isq,grad_Isq] = variance_intensity(obj,moments,sample_time)
            %  compute variance from moments
            Nsq = sum(obj.alpha.^2.*moments.*(1-moments));
            prefactor = exp(-2*obj.obs_param(3)*sample_time)*obj.obs_param(4)^2;
            Isq = prefactor*Nsq;
            if nargout == 2
                grad_Isq = prefactor*obj.alpha.^2.*(1-2*moments);
            end
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

