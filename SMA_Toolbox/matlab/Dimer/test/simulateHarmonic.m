function simulateHarmonic(x0, c, D, a, dt, N)
    % x0 - Start position (micron)
    % c - Well center (micron)
    % D - Diffusion constant (micron^2/s)
    % a - Force units of (kB*T)
    % dt - Time step (s)
    % N - number of steps
    dim =1;
    X=zeros(N,dim);
    X(1,:)=x0;
    B=sqrt(2*D); %Noise term
    A= @(x) -2*a*(x-c); % Drift term
    dt=(sqrt(2*D)/max([A(c),A(c+B),A(c-B)]))^2;
    T=N*dt;
    for n=2:N
        X(n,:) = X(n-1,:) + randn(1,1)*B*sqrt(dt) + A(X(n-1,:))*dt;
    end
    
    figure();
    plot(0,x0,'ro','DisplayName','x0: StartSite');
    hold();
    plot([0,T],[c,c],'b:','DisplayName','c: Center of well');
    plot( ((1:N)-1)*dt, X(:,1),'k.','DisplayName', 'x(t): simulated position');
    legend('location','best');
    xlabel('Time(t)');
    ylabel('Position x');
    
    S=sqrt(var(X(N/2:end,:)))
end