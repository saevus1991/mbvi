classdef BirthDeathCentral
    % Class that implements a simple birth death process
    % the central moment is used for all second order moments
    
    properties
        initial
        initial_rates
        num_reactions
        tspan
        num_moments
    end
    
    methods
        
        %% constructor
        function obj = BirthDeathCentral(initial,initial_rates,tspan)
            % store input
            obj.initial = initial;
            obj.initial_rates = initial_rates;
            obj.tspan = tspan;
            % set derived quantities
            obj.num_moments = length(initial);
            obj.num_reactions = length(initial_rates);
        end
        
        
        %% main functions
        
        % compute the moment equations
        function dydt = moments(obj,t,m,alpha_t,alpha)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
            % evaluate the derivative
            % evaluate the derivative
            dydt = zeros(2,1);
            dydt(1) = alpha(1)-alpha(2)*m(1);
            dydt(2) = alpha(1)-2*alpha(2)*m(2)+alpha(2)*m(1);
        end
        
        % compute the jacobian of the forward equation with respect to the moments
        function jacobian = forward_jacobian(obj,alpha)
            % initialize
            jacobian = zeros(2,2);
            % first row
            jacobian(1,1) = -alpha(2);
            % second row
            jacobian(2,1) = alpha(2);
            jacobian(2,2) = -2*alpha(2);
        end
        
        % compute the forward jacobian with respect to the alphas
        function jacobian = forward_jacobian_alpha(obj,m)
            % initialize
            jacobian = zeros(2,2);
            % first row
            jacobian(1,1) = 1;
            jacobian(1,2) = -m(1);
            % second row
            jacobian(2,1) = 1;
            jacobian(2,2) = m(1)-2*m(2);
        end
               
        % evalute the contribution of the scaling factors to the kl
        function alpha_eff = kl_equation(obj,alpha,c)
            % initialize
            alpha_eff = c(:,1)-alpha+alpha.*log(alpha./c(:,2));
        end
        
        % compute the propensities
        function prop = propensities(obj,m)
            % set up output
            prop = zeros(size(m,1),2);
            % birth
            prop(:,1) = 1;
            % death
            prop(:,2) = m(:,1);
        end
        
        % compute jacobian for effective alpha function
        function jacobian = prop_jacobian(obj,m)
            % initialize
            jacobian = zeros(2,2);
            % second row
            jacobian(2,1) = 1;
        end
        
        % compute the backward equation
        function deta = constraints(obj,t,eta,alpha_t,alpha,m,c)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1))';
            m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1))';
            % compute control part
            alpha_eff = obj.kl_equation(alpha,c);
            % compute Jacobian
            pj = obj.prop_jacobian(m);
            fj = obj.forward_jacobian(alpha);
            % evaluae result
            deta = pj'*alpha_eff-fj'*eta;
        end
        
        % gradient of the control
        function grad = control_gradient(obj,eta,m)
            % number of time_steps
            n_t = size(m,1);
            % initialize
            grad = zeros(n_t,2);
            % contributions lagrange multpliers
            for i = 1:n_t
                fj = obj.forward_jacobian_alpha(m(i,:)');
                grad(i,:) = eta(i,:)*fj;
            end
        end
        
    end
end

