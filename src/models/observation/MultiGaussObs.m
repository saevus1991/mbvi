classdef MultiGaussObs < handle
    % Gaussian observation model for a secnerio where several statistics
    % are observed independently with individual noises stored in a vector sigma
    
    properties
        sigma
        moment_selection
    end
    
    methods
        
        % constructor
        function obj = MultiGaussObs(sigma,num_species,observed_moments)
            obj.sigma = sigma;
            obj.moment_selection = obj.init_moment_selection(num_species,observed_moments);
        end
        
        % compute which indices have to be updated
        function moment_selection = init_moment_selection(obj,num_species,observed)
            % initialize
            moment_selection = zeros(length(observed),2);
            % linear moments
            moment_selection(:,1) = observed;
            % corresponding square moments
            moment_selection(:,2) = num_species+(observed-1).*(num_species-0.5*observed+1)+1;
        end
        
        % compute terminal condition of the backward equation
        function terminal = get_terminal(obj,moments,observation,sample_time)
            % initialize
            terminal = zeros(size(moments,1),1);
            % update observed species
            terminal(obj.moment_selection(:,1)) = (observation-moments(obj.moment_selection(:,1)))./obj.sigma.^2;
            terminal(obj.moment_selection(:,2)) = -0.5./obj.sigma.^2;
        end
               
        % compute the residual
        function residual = get_residual(obj,moments,observation,sample_time)
            % contribution of the normalizer
            residual = 0.5*length(obj.sigma)*log(2*pi)+sum(log(obj.sigma));
            % central moment contribution
            residual = residual+0.5*sum(moments(obj.moment_selection(:,2))./obj.sigma.^2);
            % observation contribution
            residual = residual+0.5*sum((observation-moments(obj.moment_selection(:,1))).^2./obj.sigma.^2);
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

