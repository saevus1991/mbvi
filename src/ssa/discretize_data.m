function [ discrete_data ] = discretize_data(t_true,x_data,t_sample,selection )
% Discretize a jump process path
%   input:
%       sample_data: 1xnum_cells cell struct containing cell arrays of the form
%           {time,reaction_count,sample_path,sample_path_integral}
%       selecion: integer array indicating which species are measured
%       delta_t: sampling time. If delta_t is a 2d vector, the second
%       component is considered as a time offset
%   output:
%       discrete_data: vector of state evaluated at sample time

%% input check
if nargin == 3
    selection = 1:size(x_data,1);
end


%% preparations 
cnt = 1;
ind = zeros(size(t_sample));

%% create data
% find indices corresponding to the sampling times
for i = 1:length(t_sample)
    while t_true(cnt+1) < t_sample(i)
        cnt = cnt + 1;
    end
    ind(i) = cnt;
end

% extract indexed data and apply noise
discrete_data = x_data(selection,ind);


end

