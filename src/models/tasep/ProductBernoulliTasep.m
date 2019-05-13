classdef ProductBernoulliTasep
    % Class that implements forward, backward and constraint functions for
    % the TASEP process wit 3 variational scaling factors
    
    properties
        num_sites
        initial
        initial_rates
        num_reactions
        tspan
        num_moments
        num_reaction
    end
    
    methods
        
        %% constructor
        function obj = ProductBernoulliTasep(num_sites,initial,initial_rates,tspan)
            obj.num_sites = num_sites;
            obj.initial = initial;
            obj.initial_rates = initial_rates;
            obj.tspan = tspan;
            obj.num_moments = length(initial);
            obj.num_reactions = length(initial_rates);
        end
        
        
        %% main functions
        
        % compute the Bernoulli closure approximation of the mean
        function dm = moments(obj,t,m,alpha_t,alpha)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
            % evalute derivative
            dm = zeros(obj.num_sites,1);
            dm(1) = alpha(1)*(1-m(1))-alpha(2)*m(1)*(1-m(2));
            dm(2:obj.num_sites-1) = alpha(2)*(m(1:obj.num_sites-2).*(1-m(2:obj.num_sites-1))-m(2:obj.num_sites-1).*(1-m(3:obj.num_sites)));
            dm(obj.num_sites) = alpha(2)*m(obj.num_sites-1)*(1-m(obj.num_sites))-alpha(3)*m(obj.num_sites);
        end
        
        % compute the jacobian of the forward equation
        function jacobian = forward_jacobian(obj,m,alpha)
            % initialize
            jacobian = zeros(obj.num_sites,obj.num_sites);
            % set diagonal
            jacobian(1,1) = -alpha(1)-alpha(2)*(1-m(2));
            for i = 2:obj.num_sites-1
                jacobian(i,i) = -alpha(2)*m(i-1)-alpha(2)*(1-m(i+1));
            end
            jacobian(end,end) = -alpha(2)*m(end-1)-alpha(3);
            % set lower diagonal
            for i = 2:obj.num_sites
                jacobian(i,i-1) = alpha(2)*(1-m(i));
            end
            % set uppder diagonal
            for i = 1:obj.num_sites-1
                jacobian(i,i+1) = alpha(2)*m(i);
            end
        end
        
        % compute jacobian with respect to alpha
        function jacobian = forward_jacobian_alpha(obj,m)
            % initialize
            jacobian = zeros(obj.num_sites,3);
            % first column
            jacobian(1,1) = (1-m(1));
            % second column
            jacobian(1,2) = -m(1)*(1-m(2));
            jacobian(2:end-1,2) = m(1:obj.num_sites-2).*(1-m(2:obj.num_sites-1))-m(2:obj.num_sites-1).*(1-m(3:obj.num_sites));
            jacobian(end,2) = m(obj.num_sites-1)*(1-m(obj.num_sites));
            % third column
            jacobian(end,3) = -m(obj.num_sites);
        end
               
        % evalute the contribution from the kl
        function alpha_eff = kl_equation(obj,alpha,c)
            % initialize
            alpha_eff = zeros(3,1);
            % stretch to correct form
            alpha_eff(1) = c(1,1)-alpha(1)+alpha(1)*log(alpha(1)/c(1,2));
            alpha_eff(2) = c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2));
            alpha_eff(end) = c(3,1)-alpha(3)+alpha(3)*log(alpha(3)/c(3,2));
        end
        
        % compute jacobian for effective alpha function
        function jacobian = prop_jacobian(obj,m)
            % initialize
            jacobian = zeros(obj.num_sites+1,obj.num_sites);
            % set diagonal
            jacobian(1,1) = -1;
            for i = 2:obj.num_sites
                jacobian(i,i) = -m(i-1);
            end
            % set sub-diagonal
            for i = 2:obj.num_sites
                jacobian(i,i-1) = (1-m(i));
            end
            jacobian(obj.num_sites+1,obj.num_sites) = 1;
            % concatenate for parameter-based partition
            jacobian = [jacobian(1,:);sum(jacobian(2:end-1,:));jacobian(end,:)];
        end
        
        % compute the backward equation
        function deta = constraints(obj,t,eta,alpha_t,alpha,m,c)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
            m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1));
            % compute control part
            alpha_eff = obj.kl_equation(alpha,c);
            % compute Jacobian
            pj = obj.prop_jacobian(m);
            fj = obj.forward_jacobian(m,alpha);
            % evaluae result
            deta = pj'*alpha_eff-fj'*eta;
        end

        % compute the propensities
        function prop = propensities(obj,m)
            prop = [1-m(:,1),m(:,1:end-1).*(1-m(:,2:end)),m(:,end)];
            % parameter partition
            prop = [prop(:,1),sum(prop(:,2:end-1),2),prop(:,end)];
        end
        
        % gradient of the control
        function grad = control_gradient(obj,eta,m)
            % number of time_steps
            n_t = size(m,1);
            % initialize
            grad = zeros(n_t,3);
            % contributions lagrange multpliers
            for i = 1:n_t
                fj = obj.forward_jacobian_alpha(m(i,:)');
                grad(i,:) = eta(i,:)*fj;
            end
        end
        
    end
end

