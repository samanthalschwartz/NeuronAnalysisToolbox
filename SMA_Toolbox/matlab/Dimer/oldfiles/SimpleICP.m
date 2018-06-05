

classdef SimpleICP
    properties
            D_A = 0.1;
            D_B = 0.1;
            D_AB = 0.05;
            dT = 0.02;
            
            SE_R = 0.020;
            SE_M = 0.020;
    end
    methods
        
        function [freeLLH, boundLLH] = computeModelLLH(obj, ObsA, ObsB)
            V_FA = 2*obj.dT*obj.D_A + obj.SE_M^2;
            V_FB = 2*obj.dT*obj.D_B + obj.SE_M^2;
%             freeLLH = -0.5*(2*log(2*pi) + log(V_FA) + log(V_FB) + diff(ObsA).^2./V_FA + diff(ObsB).^2./V_FB);
            freeLLH = log(1-normcdf(abs(diff(ObsA)),0,sqrt(V_FA))) + log(1-normcdf(abs(diff(ObsB)),0,sqrt(V_FB)));
            
            V_BA = 2*obj.dT*obj.D_AB + obj.SE_M^2;
            V_BR = obj.SE_R^2 + 2*obj.SE_M^2;
%             boundLLH = -0.5*(2*log(2*pi) + log(V_BA) + log(V_BR) + diff(ObsA).^2./V_BA + (ObsA(1:end-1)-ObsB(1:end-1)).^2./V_BR);
            boundLLH = log(1-normcdf(abs(diff(ObsA)),0,sqrt(V_BA))) + log(1-normcdf(abs(ObsA(1:end-1)-ObsB(1:end-1)),0,sqrt(V_BR)));
            
            
            disp_A = abs(ObsA(2:end)-ObsA(1:end-1));
            disp_B = abs(ObsB(2:end)-ObsB(1:end-1));
            dist = ObsA(2:end) -ObsB(2:end);
            
            SE_Dimer = sqrt(obj.SE_R^2+obj.SE_M^2);
            freeLLH = log(normpdf(abs(disp_B),0,0.5*(disp_B+disp_A)));
            boundLLH = log(normpdf(abs(dist),0,SE_Dimer));
        end
        
        function LLH = computeSequenceLLH(obj, ObsA, ObsB, state0, cps)
            [freeLLH, boundLLH] = obj.computeModelLLH(ObsA, ObsB);
            N=numel(ObsA);
            LLH=0;
            state=state0;
            ncp = numel(cps);
            cpidx=1;
            for i=1:N-1
                if cpidx<=ncp && cps(cpidx)==i
                    state = ~state;
                    cpidx = cpidx+1;
                end
                if state
                    LLH = LLH + boundLLH(i);
                else
                    LLH = LLH + freeLLH(i);
                end
            end         
        end
        
        function [H0_LLH, H1_LLH] = computeNewLLH(obj, ObsA, ObsB)
            disp_A = abs(ObsA(2:end)-ObsA(1:end-1));
            disp_B = abs(ObsB(2:end)-ObsB(1:end-1));
            dist = ObsA(2:end) -ObsB(2:end);
            
            SE_Dimer = sqrt(obj.SE_R^2+obj.SE_M^2);
            H0_LLH = log(normpdf(abs(disp_A),0,2*(disp_B+disp_A)));
            H1_LLH = log(normpdf(abs(dist),0,SE_Dimer));
            
            figure();            
            Ts=1:numel(ObsA);
            subplot(2,1,1);
            plot(Ts,ObsA,'r-');
            hold on;
            plot(Ts,ObsB,'b-');
%             for cp=cps
%                 plot([cp,cp],ylim(),'k:');
%             end
            subplot(2,1,2);
            Ts = 2:numel(ObsA);
            plot(Ts, H0_LLH, 'r-');
            hold on;
            plot(Ts, H1_LLH, 'b-');
             
        end
        
        
        function [MLE_state0, MLE_cps, MLE_LLH] = computeSequenceMLE(obj, ObsA, ObsB)
            [freeLLH, boundLLH] = obj.computeModelLLH(ObsA, ObsB);
            freeLLH_cum = [0; cumsum(freeLLH)];
            boundLLH_cum = [0; cumsum(boundLLH)];
            
            N = numel(ObsA);
            beta=log(N);
            MLE_F = zeros(N-1,1);
            MLE_B = zeros(N-1,1);
            cps_F = zeros(N-1,1);
            cps_B = zeros(N-1,1);
            for n = N-1:-1:1
               cps_F(n)=N;
               MLE_F(n) = freeLLH_cum(end) - freeLLH_cum(n);
               for m=n+1:N-1
                   temp = (freeLLH_cum(n) - freeLLH_cum(m)) + MLE_B(m) - beta;
                   if temp>MLE_F(n)
                       MLE_F(n) = temp;
                       cps_F(n) = m;
                   end
               end
               cps_B(n)=N;
               MLE_B(n) = boundLLH_cum(end) - boundLLH_cum(n);
               for m=n+1:N-1
                   temp = (boundLLH_cum(n) - boundLLH_cum(m)) + MLE_F(m) - beta;
                   if temp>MLE_B(n)
                       MLE_B(n) = temp;
                       cps_B(n) = m;
                   end
               end               
            end
            MLE_state0 = MLE_B(1) > MLE_F(1);
            if MLE_state0
                MLE_LLH = MLE_B(1);
                next_cp = cps_B(1);
            else
                MLE_LLH = MLE_F(1);               
                next_cp = cps_F(1);
            end
            state = MLE_state0;
            MLE_cps = zeros(N-1,1);
            k=1;
            
            while next_cp<N;
                MLE_cps(k) = next_cp;
                k = k + 1;
                state = ~state;
                if state
                    next_cp = cps_B(next_cp);
                else
                    next_cp = cps_F(next_cp);
                end
            end
            MLE_cps = MLE_cps(1:k-1);
        end
        
        function f = plotDimer(obj, ObsA, ObsB, cps)
            if nargin==3
                cps=[];
            end
            [freeLLH, boundLLH] = obj.computeModelLLH(ObsA, ObsB);
            
            f = figure();            
            Ts=1:numel(ObsA);
            subplot(4,1,1);
            plot(Ts,ObsA,'r-');
            hold on;
            plot(Ts,ObsB,'b-');
            for cp=cps
                plot([cp,cp],ylim(),'k:');
            end
            subplot(4,1,2);
            Ts=1:(numel(ObsA)-1);
            plot(Ts, freeLLH, 'r-');
            hold on;
            plot(Ts, boundLLH, 'b-');

            subplot(4,1,3);
            Ts=1:(numel(ObsA)-1);
            plot(Ts, cumsum(freeLLH), 'r-');
            hold on;
            plot(Ts, cumsum(boundLLH), 'b-');

            subplot(4,1,4);
            Ts=1:(numel(ObsA)-1);
            hold on;
            ratio = 1./(1+exp(cumsum(freeLLH) - cumsum(boundLLH)));
            plot(Ts, 1-ratio, 'r-');
            plot(Ts, ratio, 'b-');
            ratio = 1./(1+exp(freeLLH - boundLLH));
            plot(Ts, 1-ratio, 'r:');
            plot(Ts, ratio, 'b:');

        end
        
        function [ObsA, ObsB] = simulateDimer(obj, N, cps, state0)

            SE_DA = sqrt(2*obj.dT*obj.D_A);
            SE_DB = sqrt(2*obj.dT*obj.D_B);
            SE_DAB = sqrt(2*obj.dT*obj.D_AB);


            ObsA = zeros(N,1);
            ObsB = zeros(N,1);
            ncps = numel(cps);
            ext_cps=[1, cps, N];
            state = state0;
            for cp_idx = 1:(ncps+1)
                cp0 = ext_cps(cp_idx);
                cp1 = ext_cps(cp_idx+1);
                if state
                    SE_A = SE_DAB;
                else
                    SE_A = SE_DA;
                end
                n = cp1-cp0;
                rng = cp0+1:cp1;
                ObsA(rng) = ObsA(cp0)+cumsum(SE_A*randn(n,1))+obj.SE_M*randn(n,1) ;
                if state
                    ObsB(rng) = ObsA(rng) + sqrt(obj.SE_R^2+obj.SE_M^2)*randn(n,1);
                else
                    TrueB = ObsA(cp0)+cumsum(SE_DB*randn(n,1))+obj.SE_M*randn(n,1);
                    if cp1==N
                        ObsB(rng) = TrueB;
                    else
                        ObsB(rng) = TrueB + linspace(0,ObsA(cp1)-TrueB(end),n)';
                    end
                end

                state = ~state;
            end
           
        end 
    end %public methods

end %classdef
