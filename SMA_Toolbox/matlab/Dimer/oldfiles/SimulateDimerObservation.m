classdef SimulateDimerObservation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        D=0.1; % um/s^2
        dT=0.1; %100 fps
        dBond=0.050; %(um) Bond distance
        sigmaX=0.010; % (um)
        sigmaY=0.008; % (um)
        sigmaRegX=0.005; %(um)
        sigmaRegY=0.005; %(um)
    end
    
    methods
        function [obsTraj,trueTraj] = simulateBoundTrajectory(obj, N)
            times = obj.dT*(0:(N-1))';
            sigmaD = 2*obj.dT*obj.D;
            centerX = [0; cumsum(randn(N-1,1)*sigmaD)];
            centerY = [0; cumsum(randn(N-1,1)*sigmaD)];
            angle = 2*pi/360*randi([0,360],N,1);
            X1 = centerX + 0.5*obj.dBond*cos(angle);
            Y1 = centerY + 0.5*obj.dBond*sin(angle);
            X2 = centerX + 0.5*obj.dBond*cos(angle+pi);
            Y2 = centerY + 0.5*obj.dBond*sin(angle+pi);
            obsX1 = X1+randn(N,1)*obj.sigmaX;
            obsY1 = Y1+randn(N,1)*obj.sigmaY;
            obsX2 = X2+randn(N,1)*obj.sigmaX+rand(N,1)*obj.sigmaRegX;
            obsY2 = Y2+randn(N,1)*obj.sigmaY+rand(N,1)*obj.sigmaRegY;
            obsTraj=[times, obsX1, obsX2, obsY1,obsY2];
            trueTraj=[times, X1, X2, Y1, Y2];
            figure('Position',[10,10,500,1000]);
            subplot(3,1,1);
%             plot(obsX1,obsY1,'r-','linewidth',2);
            hold('on');
            axis('equal')
%             plot(obsX2,obsY2,'b-','linewidth',2);
            plot(X1,Y1,'r-');
            plot(X2,Y2,'b-');
            subplot(3,1,2);
            D=sqrt((X1-X2).^2+(Y1-Y2).^2);
            obsD = sqrt((obsX1-obsX2).^2+(obsY1-obsY2).^2);
            plot(times,D,'k-');
            hold('on');
            plot(times,obsD,'r-');
            subplot(3,1,3);
            H = histogram(obsD,'BinMethod','scott');
            hold('on');
            sigmaFit = 0.5*(obj.sigmaX+obj.sigmaY+obj.sigmaRegX+obj.sigmaRegY);
            plot(H.BinEdges(1:end-1)+H.BinWidth/2,sum(H.Values)*diff(normcdf(H.BinEdges,obj.dBond,sigmaFit)),'r-');
        end
        function [obsTraj,trueTraj] = simulateFreeTrajectory(obj, N)
            times = obj.dT*(0:(N-1))';
            sigmaD = 2*obj.dT*obj.D;
            X1 = [0; cumsum(randn(N-1,1)*sigmaD)];
            Y1 = [0; cumsum(randn(N-1,1)*sigmaD)];
            X2 = [0; cumsum(randn(N-1,1)*sigmaD)];
            Y2 = [0; cumsum(randn(N-1,1)*sigmaD)];
            obsX1 = X1+randn(N,1)*obj.sigmaX;
            obsY1 = Y1+randn(N,1)*obj.sigmaY;
            obsX2 = X2+randn(N,1)*obj.sigmaX+rand(N,1)*obj.sigmaRegX;
            obsY2 = Y2+randn(N,1)*obj.sigmaY+rand(N,1)*obj.sigmaRegY;
            obsTraj=[times, obsX1, obsX2, obsY1,obsY2];
            trueTraj=[times, X1, X2, Y1, Y2];
            figure('Position',[10,10,500,1000]);
            subplot(3,1,1);
%             plot(obsX1,obsY1,'r-','linewidth',2);
            hold('on');
            axis('equal')
%             plot(obsX2,obsY2,'b-','linewidth',2);
            plot(X1,Y1,'r-');
            plot(X2,Y2,'b-');
            subplot(3,1,2);
            D=sqrt((X1-X2).^2+(Y1-Y2).^2);
            obsD = sqrt((obsX1-obsX2).^2+(obsY1-obsY2).^2);
            plot(times,D,'k-');
            hold('on');
            plot(times,obsD,'r-');
            subplot(3,1,3);
            H = histogram(diff(obsD),'BinMethod','scott');
%             hold('on');
%             sigmaFit = 0.5*(obj.sigmaX+obj.sigmaY+obj.sigmaRegX+obj.sigmaRegY);
%             plot(H.BinEdges(1:end-1)+H.BinWidth/2,sum(H.Values)*diff(normcdf(H.BinEdges,obj.dBond,sigmaFit)),'r-');
        end

    end
    
end

