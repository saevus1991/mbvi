 function [times,states] = ssa_sparse(initial,tspan,generator)
 %% stochastic simulation algrotih for a Markov chain with sparse generator
 
 % initialize
 state = initial;
 time = tspan(1);
 times = [];
 states = [];
 
 % split genertor
 [rows,cols,values] = find(generator);
 
 % draw waiting times
 while time < tspan(2)
     % get column indices corresponding to the current state
     ind = (cols == state & rows ~= state);
     % select propensities
     prop = values(ind);
     total_prop = sum(prop);
     % draw waiting time
     if total_prop == 0
         times = [times,tspan(2)];
         states = [states,state];
         break;
     else
        delta_t = -log(rand)/total_prop;
     end
     % draw event
     event = randsample(1:length(prop),1,true,prop);
     % update time and state
     time = time+delta_t;
     state = rows(ind);
     state = state(event);
     % save time and state
     times = [times,time];
     states = [states,state];
 end
 
 
 end
 