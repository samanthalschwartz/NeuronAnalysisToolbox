classdef MultiLifetimeEstimator < handle
    properties
        nTracks;
        trackLengths;
        maxLength;
    end %properties
    
    methods
        function obj=MultiLifetimeEstimator(trackLengths)
            obj.trackLengths=trackLengths;
            obj.nTracks=length(obj.trackLengths);
            obj.maxLength=max(obj.trackLengths);
        end
        
        function [thetaMLE, maxLLH] = MLEsingle(obj,initial_guess)
            fprintf('Maximizing: Initial guess R:%.9f\n',initial_guess(1));
            fprintf('Initial LLH: %.9f\n',obj.singleLifetime(initial_guess));
            [thetaMLE, maxLLH, flags, stats]=fminsearch(@(theta) -obj.singleLifetime(theta), initial_guess);
            fprintf('ThetaMLE %.9f | MaxLLH: %.9f\n',thetaMLE,maxLLH);
            fprintf('Flags: %s\n',flags);
            fprintf('Stats:');
            disp(stats);
        end

        function LLH=singleLifetime(obj, theta)
            % theta=[rD]
            LLH=0;
            for track=1:obj.nTracks
                n=obj.trackLengths(track);
                LLH=LLH+ (n-1)*log(1-theta(1))+log(theta(1)); %n-1 tails and then 1 heads
            end
        end
        function plotTrackLengths(obj)
            count=histc(obj.trackLengths, 1:obj.maxLength);
            bar(count);
        end
    end

    methods (Static=true)
        function mlife_est=simulateDataset_single(n, r)
            % n - number of tracks to simulate
            % rD -  dissociation rate/timestep
            trackLengths=zeros(1,n);
            for track=1:n
                k=1;
                while( rand()>= r); k=k+1; end;
                trackLengths(track)=k;
            end
            mlife_est=MultiLifetimeEstimator(trackLengths);
        end
    end % Static Methods

end % classdef 
