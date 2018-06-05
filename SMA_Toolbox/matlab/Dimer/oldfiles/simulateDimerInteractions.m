function [ObsA, ObsB] = simulateDimerInteractions(N, cps, state0)
    D_A = 0.1;
    D_B = 0.1;
    D_AB = 0.005;
    dT = 0.002;
    
    SE_DA = sqrt(2*dT*D_A);
    SE_DB = sqrt(2*dT*D_B);
    SE_DAB = sqrt(2*dT*D_AB);

    SE_R = 0.02;
    SE_E = 0.03;
    
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
        ObsA(rng) = ObsA(cp0)+cumsum(SE_A*randn(n,1))+SE_E*randn(n,1) ;
        if state
            ObsB(rng) = ObsA(rng) + sqrt(SE_R^2+SE_E^2)*randn(n,1);
        else
            TrueB = ObsA(cp0)+cumsum(SE_DB*randn(n,1));
            if cp1==N
                ObsB(rng) = TrueB;
            else
                ObsB(rng) = TrueB + linspace(0,ObsA(cp1)-TrueB(end),n)';
            end
        end
        
        state = ~state;
    end
    figure();
    Ts=1:N;
    plot(Ts,ObsA,'r-');
    hold on;
    plot(Ts,ObsB,'b-');
    for cp=cps
        plot([cp,cp],ylim(),'k:');
    end
end

