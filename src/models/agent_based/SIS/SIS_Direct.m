classdef SIS_Direct
    % Class that implements a gene expression model with a single gene mrna and protein
    % Central moments are used for all second order moments
    
    properties
        initial
        initial_rates
        graph
        tspan
        num_nodes
        num_states
        num_reactions
    end
    
    methods
        
        %% constructor
        function obj = SIS_Direct(initial,initial_rates,graph,tspan,num_states)
            % store input
            obj.initial = initial;
            obj.initial_rates = initial_rates;
            obj.graph = graph;
            obj.tspan = tspan;
            % set derived quantities
            obj.num_nodes = length(initial);
            obj.num_states = num_states;
            obj.num_reactions = 2*obj.num_nodes;
        end
        
        
        %% main functions
        
        % compute the moment equations (this is general for agent basedmodels)
        function dydt = moments(obj,t,m,alpha_t,alpha)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1))';
            % compute meanfield infection rate
            phi_SI = alpha(1:2:end).*(1-m).*(obj.graph*m);
            % initialize
            dydt = -alpha(2:2:end).*m+phi_SI;
        end
        
        % compute the effective rate matrix as obtained from multiplication with variational scaling factors
        function jacobian = forward_jacobian(obj,m,alpha)
            % initialize
            jacobian = zeros(obj.num_nodes,obj.num_nodes);
            % filll the diagonal
            tmp = alpha(1:2:end).*(obj.graph*m);
            jacobian(1:obj.num_nodes+1:end) = -tmp-alpha(2:2:end);
            % add the rest
            jacobian = jacobian+(obj.graph.*(1-m)).*alpha(1:2:end);
        end
        
        % compute the forward jacobian with respect to the alphas
        function jacobian = forward_jacobian_alpha(obj,m)
            % set up output
            jacobian = zeros(obj.num_nodes,2*obj.num_nodes);
            % add on->off contributions
            jacobian(obj.num_nodes+1:2*obj.num_nodes+1:end) = -m;
            % add off->on contributions
            tmp = (obj.graph*m).*(1-m);
            jacobian(1:2*obj.num_nodes+1:end) = tmp;
        end
               
        % evalute the contribution of the scaling factors to the kl
        function alpha_eff = kl_equation(obj,alpha,c)
            % calculate
            alpha_eff = c(:,1)-alpha+alpha.*log(alpha./c(:,2));
            % remove nans
            alpha_eff(isnan(alpha_eff)) = 0;
        end
        
        % compute the propensities
        function prop = propensities(obj,m)
            % initialize
            prop = zeros(size(m,1),2*obj.num_nodes);
            % contributions off->on transitions
            prop(:,1:2:end) = (m*obj.graph').*(1-m);
            % contributions fromm on->off transitions
            prop(:,2:2:end) = m;
        end
        
        % compute jacobian for effective alpha function
        function jacobian = prop_jacobian(obj,m)
            % initialize
            jacobian = zeros(2*obj.num_nodes,obj.num_nodes);
            % neighbor contribution from off->on
            jacobian(1:2*obj.num_nodes+2:end) = -obj.graph*m;
            % self contribution from off ->on
            jacobian(1:2:end,:) = (1-m).*obj.graph;
            % contributions from on->off
            jacobian(2:2*obj.num_nodes+2:end) = 1;
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
            fj = obj.forward_jacobian(m,alpha);
            % evaluae result
            deta = pj'*alpha_eff-fj'*eta;
        end
        
        % gradient of the control
        function grad = control_gradient(obj,eta,m)
            % number of time_steps
            n_t = size(m,1);
            % initialize
            grad = zeros(n_t,2*obj.num_nodes);
            % contributions lagrange multpliers
            for i = 1:n_t
                % compute gradient contribution
                fj = obj.forward_jacobian_alpha(m(i,:)');
                grad(i,:) = eta(i,:)*fj;
            end
        end
        
    end
end

