classdef LognormObs < handle
    % Lognormal observation model approximated by taylor expansioon for a single observed species
    
    properties
        sigma
        moment_selection
        reset_type
    end
    
    methods
        
        % constructor
        function obj = LognormObs(sigma,num_species,observed_moment)
            obj.sigma = sigma;
            obj.moment_selection = obj.init_moment_selection(num_species,observed_moment);
            obj.reset_type = 'reduced';
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
            terminal(obj.moment_selection(1)) = (log(observation)-log(moments(obj.moment_selection(1))))/moments(obj.moment_selection(1))/obj.sigma^2;
            if strcmp(obj.reset_type,'full')
                terminal(obj.moment_selection(1)) = terminal(obj.moment_selection(1))+0.5*(obj.moment_selection(2)./obj.moment_selection(1).^3)/obj.sigma^2;
                terminal(obj.moment_selection(1)) = terminal(obj.moment_selection(1))+(log(observation)-log(moments(obj.moment_selection(1)))+1).*(obj.moment_selection(2)./obj.moment_selection(1).^3)/obj.sigma^2;
                terminal(obj.moment_selection(2)) = -0.5*(log(observation)-log(moments(obj.moment_selection(1)))+1)./obj.moment_selection(1).^2/obj.sigma^2;
            end
        end
               
        % compute the residual
        function residual = get_residual(obj,moments,observation,sample_time)
            % contribution of the normalizer
            residual = 0.5*log(2*pi)+log(obj.sigma)+log(observation);
            % mean contribution
            residual = residual+0.5*(log(observation)-log(moments(obj.moment_selection(1)))).^2/obj.sigma^2;
            % second order correction
            residual = residual+0.5*(log(observation)-log(moments(obj.moment_selection(1)))+1).*(obj.moment_selection(2)./obj.moment_selection(1).^2)/obj.sigma^2;
        end
        
        
    end
end

