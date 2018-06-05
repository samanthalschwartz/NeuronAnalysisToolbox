function  testNu(rho,D,N)
    ts=linspace(0,10,100);
    snus = arrayfun(@(t) sampleNu(t,rho,D,N), ts);
    
    nuLT = @(s) (1-2*besseli(1,rho/sqrt(D)*sqrt(s)).*besselk(1,rho/sqrt(D)*sqrt(s)))./abs(s) ;
    nus = arrayfun(@(t) laplaceInverse(nuLT,t),ts);
    
    figure();
    plot(ts,snus,'k.','DisplayName','Sampled nu');
    hold();
    plot(ts,nus,'r-','DisplayName','Computed nu');
    xlabel('Time (t)');
    ylabel('nu(t)');

end

function nu = sampleNu(t,rho,D,N)
    aR = [sqrt(rand(N,1))*rho, rand(N,1)*2*pi]; % r's are distributed as sqrt(U([0,1]))
    a= [aR(:,1).*cos(aR(:,2)), aR(:,1).*sin(aR(:,2))];
    sigma = sqrt(2*D*t);
    a = a + randn(N,2).*sigma;
    nu = sum(sum(a.^2,2)<=rho^2)/N;
end

