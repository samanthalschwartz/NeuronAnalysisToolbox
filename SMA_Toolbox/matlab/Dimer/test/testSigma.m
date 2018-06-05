function testSigma()

    muA = [50, 10];
    muB = [-3, -5];
    sigmaAx = 3;
    sigmaAy = 3;
    corr = -.9;
    varA = [sigmaAx^2 corr*sigmaAx*sigmaAy; corr*sigmaAx*sigmaAy sigmaAy^2];
    varB = [4 -.8; -.8 4];
    muC = (muA+muB)/2;
    varC = (varA+varB)/4;
    errA = chol(varA);
    errB = chol(varB);
    errC = chol(varC);

    N = 10000;
    sampleA = randn(N,2)*errA+repmat(muA, N, 1);
    sampleB = randn(N,2)*errB+repmat(muB, N, 1);
    sampleC = randn(N,2)*errC+repmat(muC, N, 1);
    trueC = .5*(sampleA+sampleB);
    figure();
    scatter(sampleA(:,1), sampleA(:,2),'r.');
    hold on;
    scatter(sampleB(:,1), sampleB(:,2),'b.');
    scatter(sampleC(:,1), sampleC(:,2),'m.');
    scatter(trueC(:,1), trueC(:,2),'g.');
    axis('equal');
    plotElipse(muA, varA, .95);
    plotElipse(muB, varB, .95);
    plotElipse(muC, varC, .95);

end

function plotElipse(mu, V, conf)
    % conf should be between 0 and 1
    deltaChiSq = chi2inv(conf,2);
    [eigVec, eigVal] = eig(V);
    C = chol(V);
    [~,order] = sort(diag(eigVal),'descend');
    axpts = zeros(2,2);
    axlen = zeros(2,1);
    for i=1:2
        axlen(i) = sqrt(eigVal(order(i),order(i)))*sqrt(deltaChiSq);
        axpts(:,i) = mu' + eigVec(:,order(i))*axlen(i);
    end
    plot([mu(1) axpts(1,1)], [mu(2), axpts(2,1)], 'k-','LineWidth',2); 
    plot([mu(1) axpts(1,2)], [mu(2), axpts(2,2)], 'k-');

    K=200; % number of points to plot in elipse
    thetas=linspace(0,2*pi,K);
    phi = mod(atan2(eigVec(2,order(1)), eigVec(1,order(1)))+2*pi,2*pi); %angle of largest eigenvector (aka major axis)
    xs = axlen(1)*cos(phi)*cos(thetas) - axlen(2)*sin(phi)*sin(thetas) + mu(1);
    ys = axlen(1)*sin(phi)*cos(thetas) + axlen(2)*cos(phi)*sin(thetas) + mu(2);
    plot(xs,ys,'k-');
    
end
