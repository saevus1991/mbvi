classdef PredatorPrey
    % Class that implements a two species predator prey model with
    % log-normal moment closure
    % Central moments are used for all second order moments
    
    properties
        initial
        initial_rates
        num_reactions
        tspan
        num_species
    end
    
    methods
        
        %% constructor
        function obj = PredatorPrey(initial,initial_rates,tspan,num_species)
            % store input
            obj.initial = initial;
            obj.initial_rates = initial_rates;
            obj.tspan = tspan;
            obj.num_species = num_species;
        end
        
        
        %% main functions
        
        % compute the moment equations
        function dydt = moments(obj,t,m,constraint_t,constraint)
            % perform interpolation
            constraint = interp_new(alpha_t(2)-alpha_t(1),constraint,t-constraint_t(1));
            % construct generator
            Q = construct_generator(obj,m,constraint);
            % evaluate derivative
            dydt = Q'*m;
        end
        
        % construct the time dependent generator
        function Q = construct_generator(obj,m,constraint)
            % create row and column indices
            tmp = repelem((2:length(m)-1)',3);
            row_index = [tmp;length(m);length(m)];
            tmp = repmat((-1:1),[length(m)-2,1])+length(m);
            col_index = [tmp;length(m)-1;length(m)];
            % create  values
            values = zeros(size(col_index));
            % contributions from birth reaction
            tmp1 = obj.initial_rates(1)*(1:length(m)-2)'.*eff_rate;
            values(3:3:end) = eff_rate;
            % contributions from death reaction
            tmp2 = obj.initial_rates(1)*(1:length(m)-2)'.*eff_rate
            for i = 1:length(m)
                col
            end
        end

        % compute the backward equation
        function deta = constraints(obj,t,eta,alpha_t,alpha,m,c)
            % perform interpolation
            alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1))';
            m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1))';
            % compute derivative
            deta = zeros(5,1);
            deta(1) = c(1,1)-alpha(1)+alpha(1)*log(alpha(1)/c(1,2))-alpha(3)*eta(2)*m(2)+(c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2)))*m(2)+(c(3,1)-alpha(3)+alpha(3)*log(alpha(3)/c(3,2)))*m(2)-eta(1)*(alpha(1)-alpha(2)*m(2))-eta(3)*(alpha(1)+alpha(2)*m(2)-2*alpha(2)*(-2*m(1)*m(2)-m(4)+(2*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/m(1)^2+(2*(m(1)*m(2)+m(4))^2)/(m(1)*m(2))-(2*(m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^3*m(2))))-eta(4)*(-alpha(3)*m(1)*m(2)+alpha(2)*m(2)^2-alpha(3)*(m(1)*m(2)+m(4))+(2*alpha(3)*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/m(1)^2+(2*alpha(3)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2))-(2*alpha(3)*(m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^3*m(2))-(2*alpha(2)*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/(m(1)*m(2))+(alpha(2)*(m(1)*m(2)+m(4))^2*(m(2)^2+m(5)))/(m(1)^2*m(2)^2))-eta(5)*(alpha(3)*m(2)+2*alpha(3)*(-m(2)^2+(2*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/(m(1)*m(2))-((m(1)*m(2)+m(4))^2*(m(2)^2+m(5)))/(m(1)^2*m(2)^2)));
            deta(2) = c(4,1)-alpha(4)+alpha(4)*log(alpha(4)/c(4,2))+alpha(2)*eta(1)*m(1)+(c(2,1)-alpha(2)+alpha(2)*log(alpha(2)/c(2,2)))*m(1)+(c(3,1)-alpha(3)+alpha(3)*log(alpha(3)/c(3,2)))*m(1)-eta(2)*(-alpha(4)+alpha(3)*m(1))-eta(3)*(alpha(2)*m(1)-2*alpha(2)*(-m(1)^2+(2*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)*m(2))-((m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)^2)))-eta(4)*(-alpha(3)*m(1)^2+alpha(2)*m(1)*m(2)+alpha(2)*(m(1)*m(2)+m(4))+(2*alpha(3)*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)*m(2))-(2*alpha(2)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2))-(alpha(3)*(m(1)^2+m(3))*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)^2)-(2*alpha(2)*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/m(2)^2+(2*alpha(2)*(m(1)*m(2)+m(4))^2*(m(2)^2+m(5)))/(m(1)*m(2)^3))-eta(5)*(alpha(4)+alpha(3)*m(1)+2*alpha(3)*(-2*m(1)*m(2)-m(4)+(2*(m(1)*m(2)+m(4))^2)/(m(1)*m(2))+(2*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/m(2)^2-(2*(m(1)*m(2)+m(4))^2*(m(2)^2+m(5)))/(m(1)*m(2)^3)));
            deta(3) = -((alpha(3)*eta(4)*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)))-eta(3)*(2*alpha(1)-(2*alpha(2)*(m(1)*m(2)+m(4))^2)/(m(1)^2*m(2)));
            deta(4) = c(2,1)+c(3,1)-alpha(2)-alpha(3)+alpha(2)*eta(1)-alpha(3)*eta(2)+alpha(2)*log(alpha(2)/c(2,2))+alpha(3)*log(alpha(3)/c(3,2))-eta(3)*(alpha(2)-2*alpha(2)*(-m(1)+(2*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)^2*m(2))))-eta(4)*(alpha(1)-alpha(4)-alpha(3)*m(1)+alpha(2)*m(2)+(2*alpha(3)*(m(1)^2+m(3))*(m(1)*m(2)+m(4)))/(m(1)^2*m(2))-(2*alpha(2)*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/(m(1)*m(2)^2))-eta(5)*(alpha(3)+2*alpha(3)*(-m(2)+(2*(m(1)*m(2)+m(4))*(m(2)^2+m(5)))/(m(1)*m(2)^2)));
            deta(5) = (alpha(2)*eta(4)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2)-eta(5)*(-2*alpha(4)+(2*alpha(3)*(m(1)*m(2)+m(4))^2)/(m(1)*m(2)^2));
        end
        
        % gradient of the control
        function grad = control_gradient(obj,eta,m)
            % evalute the contribution to the gradient
            grad = zeros(size(m,1),4);
            grad(:,1) = -eta(:,1).*m(:,1)-eta(:,3).*(m(:,1)+2.*m(:,3))-eta(:,4).*m(:,4);
            grad(:,2) = -eta(:,1).*(-m(:,1).*m(:,2)-m(:,4))-eta(:,3).*(m(:,1).*m(:,2)+m(:,4)-2.*(-m(:,1).*(m(:,1).*m(:,2)+m(:,4))+((m(:,1).^2+m(:,3)).*(m(:,1).*m(:,2)+m(:,4)).^2)./(m(:,1).^2.*m(:,2))))-eta(:,4).*(m(:,2).*(m(:,1).*m(:,2)+m(:,4))-((m(:,1).*m(:,2)+m(:,4)).^2.*(m(:,2).^2+m(:,5)))./(m(:,1).*m(:,2).^2));
            grad(:,3) = -eta(:,2).*(m(:,1).*m(:,2)+m(:,4))-eta(:,4).*(-m(:,1).*(m(:,1).*m(:,2)+m(:,4))+((m(:,1).^2+m(:,3)).*(m(:,1).*m(:,2)+m(:,4)).^2)./(m(:,1).^2.*m(:,2)))-eta(:,5).*(m(:,1).*m(:,2)+m(:,4)+2.*(-m(:,2).*(m(:,1).*m(:,2)+m(:,4))+((m(:,1).*m(:,2)+m(:,4)).^2.*(m(:,2).^2+m(:,5)))./(m(:,1).*m(:,2).^2)));
            grad(:,4) = eta(:,2).*m(:,2)+eta(:,4).*m(:,4)-eta(:,5).*(m(:,2)-2.*m(:,5));
            % adapt sign
            grad = -grad;
        end
        
    end
end

