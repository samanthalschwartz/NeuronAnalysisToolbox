function testGaussianAnglePDF(mu, sigma)
    N=1e4;
    ps = randn(N,2)*chol(sigma)+mu;
    angles = mod(atan2(ps(:,2),ps(:,1)),2*pi);
    thetas=linspace(0,2*pi,500);
    pdf_thetas = plotGaussianAnglePDF(thetas, mu, sigma);
    [V,E]=eig(sigma,'vector');
    [E,s] = sort(E,1,'descend');
    V=V(:,s);
    figure();
    subplot(1,2,1);
    plot(ps(:,1),ps(:,2),'k.');
    hold();
    plot(mu(1),mu(2),'rs','MarkerFaceColor','r','DisplayName','Center mu');
    plot(0,0,'gd','MarkerFaceColor','g','DisplayName','Origin');
    plot(mu(1)+[0,3*sqrt(E(1))*V(1,1)], mu(2)+[0,3*sqrt(E(1))*V(2,1)], '-','LineWidth',2,'Color',[.8,.6,0],'DisplayName','V(1)');
    plot(mu(1)+[0,3*sqrt(E(2))*V(1,2)], mu(2)+[0,3*sqrt(E(2))*V(2,2)], '-','LineWidth',1,'Color',[.8,.6,0],'DisplayName','V(2)');
    axis('equal');
    subplot(1,2,2);
    hold(); 
    plot(thetas,pdf_thetas,'r-','DisplayName','Computed angle pdf');
    xlim([0,2*pi]);
    H=histogram(angles,'Normalization','pdf');
    H.DisplayName='Sampled angle pdf';
    ylabel('Probability');
    legend('Location','Best');

end

function p = gaussianAnglePDF(theta, mu, sigma)
    % Mardia 1972 (p.52)
    % Covarance matrix is [sigma1^2 sigma1*sigma2*rho; sigma1*sigma2*rho sigma2^2]
    sigmaX = sqrt(sigma(1,1));
    sigmaY = sqrt(sigma(2,2));
    rho = sigma(1,2)/(sigmaX*sigmaY);

    sqrtdetsigma = sigmaX*sigmaY*sqrt(1-rho^2); % det(sigma)^{-1/2}
    detsigma = sqrtdetsigma^2;
    costheta = cos(theta);
    sintheta = sin(theta);

    C = (sigma(2,2)*costheta.^2 - sigma(2,1)*sin(2*theta) + sigma(1,1)*sintheta.^2)/detsigma;
    sqrtC = sqrt(C);
    D = (1/detsigma)*(mu(1)*(sigma(2,2)*costheta - sigma(2,1)*sintheta) + mu(2)*(sigma(1,1)*sintheta - sigma(1,2)*costheta))./sqrtC;
    G = (1/sqrtdetsigma)*(mu(1)*sintheta-mu(2)*costheta)./sqrtC;
    H = sigma(2,2)*mu(1)^2 - 2*sigma(1,2)*mu(1)*mu(2) + sigma(1,1)*mu(2)^2;
    pdfmu = (1/(2*pi*sqrtdetsigma)) * exp(-(1/(2*detsigma))*H);
    cdfD = .5*(1+erf(1/sqrt(2)*D));
    pdfG = 1/sqrt(2*pi)*exp(-.5*G.^2);
    Q = pdfmu+(1/sqrtdetsigma)*D.*cdfD.*pdfG;
    p = Q./C;
end

function p = plotGaussianAnglePDF(theta, mu,sigma)
    % Mardia 1972 (p.52)
    % Covarance matrix is [sigma1^2 sigma1*sigma2*rho; sigma1*sigma2*rho sigma2^2]
    N=numel(theta);
    sigmaX = sqrt(sigma(1,1));
    sigmaY = sqrt(sigma(2,2));
    rho = sigma(1,2)/(sigmaX*sigmaY);

    sqrtdetsigma = sigmaX*sigmaY*sqrt(1-rho^2); % det(sigma)^{-1/2}
    detsigma = sqrtdetsigma^2;
    costheta = cos(theta);
    sintheta = sin(theta);

    C = (sigma(2,2)*costheta.^2 - sigma(2,1)*sin(2*theta) + sigma(1,1)*sintheta.^2)/detsigma;
    sqrtC = sqrt(C);
    D = (1/detsigma)*(mu(1)*(sigma(2,2)*costheta - sigma(2,1)*sintheta) + mu(2)*(sigma(1,1)*sintheta - sigma(1,2)*costheta))./sqrtC;
    G = (1/sqrtdetsigma)*(mu(1)*sintheta-mu(2)*costheta)./sqrtC;
    H = sigma(2,2)*mu(1)^2 - 2*sigma(1,2)*mu(1)*mu(2) + sigma(1,1)*mu(2)^2;
%     pdfmu = mvnpdf([0,0],mu,sigma);
    pdfmu = (1/(2*pi*sqrtdetsigma)) * exp(-(1/(2*detsigma))*H);
%     cdfD = normcdf(D,0,1);
    cdfD = .5*(1+erf(1/sqrt(2)*D));
%     pdfG = normpdf(G,0,1);
    pdfG = 1/sqrt(2*pi)*exp(-.5*G.^2);
    Q = pdfmu+(1/sqrtdetsigma)*D.*cdfD.*pdfG;
    p = Q./C;
    [V,E]=eig(sigma);
    zetas=mod(atan2(V(2,:),V(1,:)),2*pi);
    phi_mu = mod(atan2(mu(2),mu(1)),2*pi);
    phi_mu_p = mod(atan2(mu(2),mu(1))+pi,2*pi);
    d_mu = sqrt(sum(mu.^2));
    m1=1/(d_mu+E(1,1)) *(mod(-phi_mu,2*pi)*d_mu +zetas(1)*E(1,1));
    m2=1/(d_mu+E(1,1)) *(mod(-phi_mu,2*pi)*d_mu + mod(zetas(1)-pi,2*pi)*E(1,1));
    w1 = mu+V(:,1)
    alpha_w1 = mod(atan2(w1(2),w1(1)),2*pi);

    [~,maxidx] = max(p);
    [~,minidx] = min(p);
    maxtheta = theta(maxidx);
    mintheta = theta(minidx);
    max1bounds = [maxtheta-pi/4, maxtheta+pi/4];
    min1bounds = [mintheta-pi/4, mintheta+pi/4];
    max2bounds = [maxtheta+pi-pi/4, maxtheta+pi+pi/4];
    min2bounds = [mintheta+pi-pi/4, mintheta+pi+pi/4];
  
    max1start = theta(maxidx);
    max2start = theta(mod(maxidx+N/2,N));
    min1start = theta(minidx);
    min2start = theta(mod(minidx+N/2,N));

    [max1,max1val] = fmincon(@(theta) -gaussianAnglePDF(theta,mu,sigma), max1start, [], [],[],[],max1bounds(1),max1bounds(2));
    [max2,max2val] = fmincon(@(theta) -gaussianAnglePDF(theta,mu,sigma), max2start, [], [],[],[],max2bounds(1),max2bounds(2));
    [min1,min1val] = fmincon(@(theta) gaussianAnglePDF(theta,mu,sigma), min1start, [], [],[],[],min1bounds(1),min1bounds(2));
    [min2,min2val] = fmincon(@(theta) gaussianAnglePDF(theta,mu,sigma), min2start, [], [],[],[],min2bounds(1),min2bounds(2));
    max1val = -max1val;
    max2val = -max2val;
    max2 = mod(max2,2*pi);
    min2 = mod(min2,2*pi);
    figure('Position',[10,10,800,600]);
    hold();
%     plot(theta,C,'r-','DisplayName','C(theta)');
%     plot(theta,D,'g--','DisplayName','D(theta)');
%     plot(theta,G,'b--','DisplayName','G(theta)');
%     plot(theta,cdfD,'g-','DisplayName','F_D(theta)');
%     plot(theta,pdfG,'b-','DisplayName','f_G(theta)');
%     plot(theta,Q,'k:','DisplayName','Q(theta)');
    plot(theta,p,'k-','DisplayName','p(theta)');
    yl=ylim();
    plot(zetas([1,1]),yl,'-','Color',[.8,.6,0],'DisplayName','phi(zeta1)');
    plot(zetas([2,2]),yl,':','Color',[.8,.6,0],'DisplayName','phi(zeta2)');
     plot(phi_mu([1,1]),yl,'m-','DisplayName','phi(mu)');
     plot(phi_mu_p([1,1]),yl,'m--','DisplayName','-phi(mu)');
     plot(alpha_w1([1,1]),yl,'b-','DisplayName','alpha(w1)');
     plot(max1([1,1]),yl,'-','Color',[.8,.4,.0],'DisplayName',sprintf('max1:[%.3f rad | %.4g deg] %.3g',max1,max1*180/pi,max1val));
     plot(max2([1,1]),yl,'--','Color',[.8,.4,.0],'DisplayName',sprintf('max2:[%.3f rad | %.4g deg] %.3g',max2,max2*180/pi,max2val));
     plot(min1([1,1]),yl,'-','Color',[0,.4,.8],'DisplayName',sprintf('min1:[%.3f rad | %.4g deg] %.3g',min1,min1*180/pi,min1val));
     plot(min2([1,1]),yl,'--','Color',[0,.4,.8],'DisplayName',sprintf('min2:[%.3f rad | %.4g deg] %.3g',min2,min2*180/pi,min2val));
%     plot(m1([1,1]),yl,'-','Color',[.8,.1,.8],'DisplayName','m1');
%     plot(m2([1,1]),yl,'--','Color',[.8,.1,.8],'DisplayName','m1');
    legend('location','best');
    xlim([0,2*pi]);
end

