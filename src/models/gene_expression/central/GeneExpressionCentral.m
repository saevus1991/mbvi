classdef GeneExpressionCentral
    % Class that implements a gene expression model with a single gene mrna and protein
    % Central moments are used for all second order moments
    
    properties
        initial
        initial_rates
        num_reactions
        tspan
        num_moments
    end
    
    methods
        
        %% constructor
        function obj = GeneExpressionCentral(initial,initial_rates,tspan)
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
            dydt = zeros(9,1);
            dydt(1) = alpha(1)*(1-m(1))-alpha(2)*m(1);
            dydt(2) = alpha(3)*m(1)-alpha(4)*m(2);
            dydt(3) = alpha(5)*m(2)-alpha(6)*m(3);
            dydt(4) = -alpha(1)*m(1)+alpha(1)-2*alpha(1)*m(4)-2*alpha(2)*m(4)+alpha(2)*m(1);
            dydt(5) = -alpha(1)*m(5)-alpha(2)*m(5)+alpha(3)*m(4)-alpha(4)*m(5);
            dydt(6) = -alpha(1)*m(6)-alpha(2)*m(6)+alpha(5)*m(5)-alpha(6)*m(6);
            dydt(7) = 2*alpha(3)*m(5)+alpha(3)*m(1)-2*alpha(4)*m(7)+alpha(4)*m(2);
            dydt(8) = alpha(3)*m(6)-alpha(4)*m(8)+alpha(5)*m(7)-alpha(6)*m(8);
            dydt(9) = 2*alpha(5)*m(8)+alpha(5)*m(2)-2*alpha(6)*m(9)+alpha(6)*m(3);
        end
        
        % compute the jacobian of the forward equation with respect to the moments
        function jacobian = forward_jacobian(obj,alpha)
            % initialize
            jacobian = zeros(9,9);
            % first row
            jacobian(1,1) = -alpha(1)-alpha(2);
            % second row
            jacobian(2,1) = alpha(3);
            jacobian(2,2) = -alpha(4);
            % third row
            jacobian(3,2) = alpha(5);
            jacobian(3,3) = -alpha(6);
            % fourth row
            jacobian(4,1) = -alpha(1)+alpha(2);
            jacobian(4,4) = -2*alpha(1)-2*alpha(2);
            % fifth row
            jacobian(5,4) =  alpha(3);
            jacobian(5,5) =  -alpha(1)-alpha(2)-alpha(4);
            % sixth row 
            jacobian(6,5) =  alpha(5);
            jacobian(6,6) =  -alpha(1)-alpha(2)-alpha(6);
            % seventh row
            jacobian(7,1) = alpha(3);
            jacobian(7,2) = alpha(4);
            jacobian(7,5) = 2*alpha(3);
            jacobian(7,7) = -2*alpha(4);
            % eighth row
            jacobian(8,6) =  alpha(3);
            jacobian(8,7) =  alpha(5);
            jacobian(8,8) =  -alpha(4)-alpha(6);
            % nineth row
            jacobian(9,2) = alpha(5);
            jacobian(9,3) = alpha(6);
            jacobian(9,8) = 2*alpha(5);
            jacobian(9,9) = -2*alpha(6);
        end
        
        % compute the forward jacobian with respect to the alphas
        function jacobian = forward_jacobian_alpha(obj,m)
            % initialize
            jacobian = zeros(9,6);
            % first row
            jacobian(1,1) = (1-m(1));
            jacobian(1,2) = -m(1);
            % second row
            jacobian(2,3) = m(1);
            jacobian(2,4) = -m(2);
            % third row
            jacobian(3,5) = m(2);
            jacobian(3,6) = -m(3);
            % fourth row
            jacobian(4,1) = -m(1)+1-2*m(4);
            jacobian(4,2) = -2*m(4)+m(1);
            % fifth row
            jacobian(5,1) = -m(5);
            jacobian(5,2) =  -m(5);
            jacobian(5,3) =  m(4);
            jacobian(5,4) = -m(5);
            % sixth row 
            jacobian(6,1) = -m(6);
            jacobian(6,2) = -m(6);
            jacobian(6,5) = m(5);
            jacobian(6,6) = -m(6);
            % seventh row
            jacobian(7,3) = 2*m(5)+m(1);
            jacobian(7,4) =-2*m(7)+m(2);
            % eighth row
            jacobian(8,3) = m(6);
            jacobian(8,4) = -m(8);
            jacobian(8,5) = m(7);
            jacobian(8,6) = -m(8);
            % nineth row
            jacobian(9,5) = 2*m(8)+m(2);
            jacobian(9,6) = -2*m(9)+m(3);
        end
               
        % evalute the contribution of the scaling factors to the kl
        function alpha_eff = kl_equation(obj,alpha,c)
            % initialize
            alpha_eff = c(:,1)-alpha+alpha.*log(alpha./c(:,2));
        end
        
        % compute the propensities
        function prop = propensities(obj,m)
            % set up output
            prop = zeros(size(m,1),6);
            % gene on
            prop(:,1) = (1-m(:,1));
            % gene off
            prop(:,2) = m(:,1);
            % translation
            prop(:,3) = m(:,1);
            % mrna degradation
            prop(:,4) = m(:,2);
            % transcription
            prop(:,5) = m(:,2);
            % protein degradation
            prop(:,6) = m(:,3);
        end
        
        % compute jacobian for effective alpha function
        function jacobian = prop_jacobian(obj,m)
            % initialize
            jacobian = zeros(6,9);
            % first row
            jacobian(1,1) = -1;
            % second row
            jacobian(2,1) = 1;
            % third row
            jacobian(3,1) = 1;
            % fourth row
            jacobian(4,2) = 1;
            % fifth row
            jacobian(5,2) = 1;
            % sixth row
            jacobian(6,3) = 1;
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
            grad = zeros(n_t,6);
            % contributions lagrange multpliers
            for i = 1:n_t
                fj = obj.forward_jacobian_alpha(m(i,:)');
                grad(i,:) = eta(i,:)*fj;
            end
        end
        
    end
end

