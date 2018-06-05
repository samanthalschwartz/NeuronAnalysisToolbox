function   testMu(r0, rho,D,N)
    ts=logspace(-3,2,1000);
    smus = arrayfun(@(t) sampleMu(t,r0,rho,D,N), ts);
    mus = arrayfun(@(t) laplaceInverse(@(s) computeLTMu(s,r0,rho,D),t),ts);
    
    figure();
    plot(ts,smus,'k.','DisplayName','Sampled mu');
    hold();
    plot(ts,mus,'r-','DisplayName','Computed mu');
    ax=gca();
    ax.XScale='log';
    xlabel('Time (t)');
    ylabel('mu(t)');
end
function ltmu=computeLTMu(s,r0,rho,D)
    omega = sqrt(s/D);
    if r0<rho % start inside trap
        ltmu = (1-omega*rho*besseli(0,omega*r0)*besselk(1,omega*rho))/s;
    else
        ltmu = (omega*rho*besseli(1,omega*rho)*besselk(0,omega*r0))/s;
    end
end
function mu = sampleMu(t,r0,rho,D,N)
    theta = rand(N,1)*2*pi; % r's are distributed as sqrt(U([0,1]))
    a= [r0.*cos(theta), r0.*sin(theta)];
    sigma = sqrt(2*D*t);
    a = a + randn(N,2).*sigma;
    mu = sum(sum(a.^2,2)<=rho^2)/N;
end