classdef GaussObs < handle
    % Gaussian observation model for a single observed species
    
    properties
        sigma
        moment_selection
    end
    
    methods
        
        % constructor
        function obj = GaussObs(sigma,num_species,observed_moment)
            obj.sigma = sigma;
            obj.moment_selection = obj.init_moment_selection(num_species,observed_moment);
        end
        
        % compute which indices have to be updated
        function moment_selection = init_moment_selection(obj,num_species,observed) 
            ind1 = observed;
            ind2 = num_species+(observed-1).*(num_species-0.5*observed+1)+1;
            moment_selection = [ind1,ind2];
        end
        
        % compute terminal condition of the backward equation
        function terminal = get_terminal(obj,moments,observation,sample_time)
            % initialize
            terminal = zeros(size(moments,1),1);
            % update observed species
            terminal(obj.moment_selection(1)) = (observation-moments(obj.moment_selection(1)))/obj.sigma^2;
            terminal(obj.moment_selection(2)) = -0.5/obj.sigma^2;
        end
               
        % compute the residual
        function residual = get_residual(obj,moments,observation,sample_time)
            % contribution of the normalizer
            residual = 0.5*log(2*pi)+log(obj.sigma);
            % central moment contribution
            residual = residual+0.5/obj.sigma^2*moments(obj.moment_selection(2));
            % observation contribution
            residual = residual+0.5/obj.sigma^2*(observation-moments(obj.moment_selection(1))).^2;
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

