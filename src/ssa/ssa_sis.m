 function [times,states] = ssa_sis(system)
 %% stochastic simulation algrotih for an sis disease model on a graph
 
 % initialize
 state = system.initial;
 time = system.t_min;
 times = [];
 states = [];
 
 % draw waiting times
 while time < system.t_max
     % calculate number of active neighbors per state
     neighbors = system.graph*state;
     % calculate propensity
     prop = system.rates(1)*(1-state).*neighbors+state*system.rates(2);
     total_prop = sum(prop);
     % check if disease has died out
     if total_prop == 0
         times = [times,system.t_max];
         states = [states,state(:,end)];
         break;
     end
     % draw waiting time
     delta_t = -log(rand)/total_prop;
     % draw event
     event = randsample(1:length(prop),1,true,prop);
     % update time and state
     time = time+delta_t;
     state(event) = rem(state(event)+1,2);
     % save time and state
     times = [times,time];
     states = [states,state];
 end
 
 
 end
 