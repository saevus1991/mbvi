function [ measurement ] = simulate_measurement(t_true,x_data,noise,t_sample,selection )
% Extract discrete, noisy measurements from complete sample data
%   input:
%       sample_data: 1xnum_cells cell struct containing cell arrays of the form
%           {time,reaction_count,sample_path,sample_path_integral}
%       selecion: integer array indicating which species are measured
%       noise: function handle describing the additive noise
%       delta_t: sampling time. If delta_t is a 2d vector, the second
%       component is considered as a time offset
%   output:
%       measurement: 2x1 cell containing the discrete time vector and a
%       vector of measurements

%% input check
if nargin == 4
    selection = 1:size(x_data,1);
end

% if length(delta_t) == 2
%     offset = delta_t(2);
%     delta_t = delta_t(1);
% else
%     offset = 0;
% end

%% preparations 
%t_sample = (min(t_true)+offset):delta_t:max(t_true);
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
measurement = x_data(selection,ind);
measurement = noise(measurement);


end

