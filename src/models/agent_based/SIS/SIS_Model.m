classdef SIS_Model
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
        function obj = SIS_Model(initial,initial_rates,graph,tspan,num_states)
            % store input
            obj.initial = initial;
            obj.initial_rates = initial_rates;
            obj.graph = graph;
            obj.tspan = tspan;
            % set derived quantities
            obj.num_nodes = length(initial)/num_states;
            obj.num_states = num_states;
            obj.num_reactions = obj.num_nodes*obj.num_states^2;
        end
        
        
        %% main functions
        
        % compute the moment equations (this is general for agent basedmodels)
        function dydt = moments(obj,t,m,alpha_t,alpha)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));
            % set up derivative vector
            dydt = zeros(obj.num_states*obj.num_nodes,1);
            % compute effective rates
            local_rates = obj.compute_local_rates(m);
            eff_rates = obj.compute_effective_rates(local_rates,alpha);
            % evaluate derivative node-wise
            for i = 1:obj.num_nodes
                % compute local derivative and store
                local_state = m((i-1)*obj.num_states+1:i*obj.num_states);
                dydt((i-1)*obj.num_states+1:i*obj.num_states) = eff_rates{i}'*local_state;
            end
        end
        
        % compute the effective rate matrix for node i (here specific for SIS )
        function local_rates = compute_local_rates(obj,m)
            % set up 
            local_rates = cell(obj.num_nodes,1);
            % iterate over nodes
            for i = 1:obj.num_nodes
                % get parents of the node
                parents = find(obj.graph(:,i));
                % compute sumary statistic
                ind = (parents-1)*obj.num_states+2;
                infected_neighbors = sum(m(ind));
                % construct Q
                Q_loc = [0,infected_neighbors;1,0];
                % store 
                local_rates{i} = Q_loc;
            end
        end
        
        % compute thelocal rate matrix for node i ( for general models )
        function eff_rates = compute_effective_rates(obj,local_rates,alpha)
            % set up 
            eff_rates = cell(obj.num_nodes,1);
            % iterate over nodes
            for i = 1:obj.num_nodes
                % construct Q
                Q_loc = local_rates{i};
                % compute local effective rate matrix
                alpha_loc = alpha((i-1)*obj.num_states^2+1:i*obj.num_states^2);
                Q_eff = Q_loc.*reshape(alpha_loc,[obj.num_states,obj.num_states])';
                Q_eff(1:obj.num_states+1:end) = 0;
                exit_rates = sum(Q_eff,2);
                Q_eff(1:obj.num_states+1:end) = -exit_rates;
                % store
                eff_rates{i} = Q_eff;
            end
        end
        
        % compute the effective rate matrix as obtained from multiplication with variational scaling factors
        function jacobian = forward_jacobian(obj,m,alpha,eff_rates)
            % initialize
            jacobian = zeros(obj.num_states*obj.num_nodes,obj.num_states*obj.num_nodes);
            % compute contributions from same node
            for i = 1:obj.num_nodes
                jacobian((i-1)*obj.num_states+1:i*obj.num_states,(i-1)*obj.num_states+1:i*obj.num_states) = eff_rates{i}';
            end
            % compute contributions of children
            for i = 1:obj.num_nodes
                parents = find(obj.graph(:,i));
                ind = (parents-1)*obj.num_states+2;
                alpha_loc = alpha((i-1)*obj.num_states^2+1:i*obj.num_states^2);
                alpha_loc = reshape(alpha_loc,[obj.num_states,obj.num_states]);
                local_state = m((i-1)*obj.num_states+1:i*obj.num_states);
                for j = parents'
                    % compute contribution for j in of state
                    infected_neighbors = sum(m(ind))-m((j-1)*obj.num_states+1);
                    Q_loc = [0,infected_neighbors;1,0];
                    Q_eff = Q_loc.*alpha_loc;
                    Q_eff(1,1) = -Q_eff(1,2);
                    Q_eff(2,2) = -Q_eff(2,1);
                    jacobian((i-1)*obj.num_states+1:i*obj.num_states,(j-1)*obj.num_states+1) = Q_eff'*local_state;
                    % compute contributon for j in on state
                    infected_neighbors = infected_neighbors+1;
                    Q_loc = [0,infected_neighbors;1,0];
                    Q_eff = Q_loc.*alpha_loc;
                    Q_eff(1,1) = -Q_eff(1,2);
                    Q_eff(2,2) = -Q_eff(2,1);
                    jacobian((i-1)*obj.num_states+1:i*obj.num_states,(j-1)*obj.num_states+2) = Q_eff'*local_state;                
                end
            end
        end
        
        % compute the forward jacobian with respect to the alphas
        function jacobian = forward_jacobian_alpha(obj,m,local_rates)
            % prepare sparse construction
            row_index = [];
            col_index = [];
            values = [];
            % iterate over columns
            for i = 1:obj.num_nodes
                % get the local state
                local_state = m((i-1)*obj.num_states+1:i*obj.num_states);
                for j = 1:obj.num_states
                    for k = 1:obj.num_states
                        if k ~= j
                            % update column
                            col_index_incr = ones(2,1)*(i-1)*obj.num_states^2+(j-1)*obj.num_states+k;
                            col_index = [col_index;col_index_incr];
                            % update row
                            row_index_incr = (i-1)*obj.num_states+[min(j,k);max(j,k)];
                            row_index = [row_index;row_index_incr];
                            % update values
                            values_incr = sign([j-k;k-j])*local_rates{i}(k,j)*local_state(k);
                            values = [values;values_incr];
                        end
                    end
                end
            end
            % construct the jacobian as sparse matrix
            jacobian = sparse(row_index,col_index,values,obj.num_nodes*obj.num_states,obj.num_nodes*obj.num_states.^2);
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
            % number of time steps
            n_t = size(m,1);
            % set up output
            prop = zeros(n_t,obj.num_nodes*obj.num_states^2);
            % iterate over time steps
            for i = 1:n_t
                % get current m
                m_tmp = m(i,:)';
                % update local reates
                local_rates = obj.compute_local_rates(m_tmp);
                % construct vector of local rates
                rates = zeros(obj.num_nodes*obj.num_states^2,1);
                for j = 1:obj.num_nodes
                    Q_loc = local_rates{j}';
                    rates((j-1)*obj.num_states^2+1:j*obj.num_states^2) = Q_loc(:);
                end
                prop(i,:) = rates.*repelem(m_tmp,obj.num_states);
            end
        end
        
        % compute jacobian for effective alpha function
        function jacobian = prop_jacobian(obj,m,local_rates)
            % initialize
            row_index = [];
            col_index = [];
            values = [];
            % iterate over rows
            for i = 1:obj.num_nodes
                for j = 1:obj.num_states
                    for k = 1:obj.num_states
                        if k ~= j
                            % compute current row index
                            row_incr = (i-1)*obj.num_states^2+(j-1)*obj.num_states+k;
                            % column element from the current node
                            col_incr = (i-1)*obj.num_states+j;
                            val_incr = local_rates{i}(j,k);
                            % check for parents contributions
                            parents = find(obj.graph(:,i));
                            ind = (parents-1)*obj.num_states+2;
                            for l = parents'
                                % compute contribution for j in of state
                                infected_neighbors = sum(m(ind))-m((j-1)*obj.num_states+1);
                                Q_loc_1 = [0,infected_neighbors;1,0];
                                % compute contributon for j in on state
                                infected_neighbors = infected_neighbors+1;
                                Q_loc_2 = [0,infected_neighbors;1,0];
                                % increment 
                                col_incr = [col_incr;(l-1)*obj.num_states+[1;2]];
                                val_incr = [val_incr;m((i-1)*obj.num_states+j)*[Q_loc_1(j,k);Q_loc_2(j,k)]];
                            end
                            row_index = [row_index;row_incr*ones(size(val_incr))];
                            col_index = [col_index;col_incr];
                            values = [values;val_incr];
                        end
                    end
                end
            end
            % construct sparse jacobian
            jacobian = sparse(row_index,col_index,values,obj.num_nodes*obj.num_states^2,obj.num_nodes*obj.num_states);
        end
        
        % compute the backward equation
        function deta = constraints(obj,t,eta,alpha_t,alpha,m,c)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1))';
            m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1))';
            % compute control part
            alpha_eff = obj.kl_equation(alpha,c);
            % compute local rates 
            local_rates = obj.compute_local_rates(m);
            eff_rates = obj.compute_effective_rates(local_rates,alpha);
            % compute Jacobian
            pj = obj.prop_jacobian(m,local_rates);
            fj = obj.forward_jacobian(m,alpha,eff_rates);
            % evaluae result
            deta = pj'*alpha_eff-fj'*eta;
        end
        
        % gradient of the control
        function grad = control_gradient(obj,eta,m)
            % number of time_steps
            n_t = size(m,1);
            % initialize
            grad = zeros(n_t,obj.num_nodes*obj.num_states^2);
            % contributions lagrange multpliers
            for i = 1:n_t
                % get current m
                m_tmp = m(i,:)';
                % update local rates
                local_rates = obj.compute_local_rates(m_tmp);
                % compute gradient contribution
                fj = obj.forward_jacobian_alpha(m_tmp,local_rates);
                grad(i,:) = eta(i,:)*fj;
            end
        end
        
    end
end

