function [f_t] =interp_new(delta_t,f,t)
% Compute a linear interpolation of the function f that is given on a time
% grid with constant step delta_t

ind = floor(t/delta_t)+1;
t_lower = delta_t*(ind-1);
f_t = f(ind,:)+(t-t_lower)/delta_t*(f(ind+1,:)-f(ind,:));

end

