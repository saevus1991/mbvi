classdef Bootstrap_Filter
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_particles
        time_s
        simulator               % a simulation angine that can generate trajectories from a given initial distribution
        observations            % values of observation
        observation_times       % liklihood model
        filter_norm             %
        trajectories
    end
    
    methods
        
        %% constructor and setup section
        
        % constructor method
        function obj = Bootstrap_Filter(num_particles)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.simulator = simulator;
            obj.num_intervals = length(sample_times)+1;
        end
        
        %% simulation section
        
        % perform the forward filter
        function forward_filter(obj)
            % 
        end
        
        % bootstrap filter update
        function outputArg = filter_step(obj)
            %
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

