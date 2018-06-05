function  [survival_prob, times] = simulateDiskCapture(r, D, T, rho, lambda, Ntrials, simT)
    % r - initial distance (um)
    % D - diffusion constant (um^2/s)
    % T - total duration (s)
    % rho - capture radius (um)
    % lambda - caputre rate when within capture disk (s^-1)
    % Ntrials - number of trials to simulate
    % simT  - simulation time step (s) [default = T/1000]
    if nargin<7
        Nsteps = 1000;
    else
        Nsteps = ceil(T/simT);
    end
    simT = T/(Nsteps-1);
    survivors = true(Ntrials, Nsteps);
    pos = zeros(Ntrials,2 ,Nsteps);
    %initialize positions
    thetas = rand(Ntrials,1)*2*pi;
    pos(:,1,1) = cos(thetas)*r;
    pos(:,2,1) = sin(thetas)*r;
    stepsize = sqrt(2*D*simT);
    capture_prob = 1-exp(-lambda*simT);
    for n=2:Nsteps
        pos(:,:,n) = pos(:,:,n-1) + randn(Ntrials,2)*stepsize;
        dists = sqrt(sum(pos(:,:,n).^2,2));
        survivors(:,n) = survivors(:,n-1) & (dists>rho | rand(Ntrials,1)>capture_prob);
    end
    
    times = linspace(0, T, Nsteps);
    survival_prob = mean(survivors,1);
    figure('Position', [10,10,600, 900]);
    subplot(2,1,1);
    plot(times, survival_prob, 'r-');
    xlabel('Time (s)');
    ylabel('Survival probability');
    ylim([0,1]);
    
    subplot(2,1,2);
    plot(pos(:,1,end), pos(:,2,end),'r.', 'MarkerSize', 2, 'DisplayName', 'Final position');
    hold();
    plot(pos(:,1,1), pos(:,2,1),'k.', 'MarkerSize', 2, 'DisplayName', 'Initial position');
    phis=linspace(0,2*pi,100);
    plot(cos(phis)*rho, sin(phis)*rho, 'b-', 'DisplayName', 'Capture radius (rho)');
    axis('equal');
    legend('Location','best')
    xlabel('pos (um)');
    ylabel('pos (um)');
end

