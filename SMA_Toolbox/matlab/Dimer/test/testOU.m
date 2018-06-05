
mu=10;
sigma = 3;
gamma=10;
x0=-10;
T=logspace(-6,2,1000);
ou = OUProcess(mu,sigma,gamma,x0);
XS = ou.sampleGillespie(T,10);
llh= ou.sampleLLH(T,XS,x0,0);
[M,V]= ou.sampleGaussianDistribution(T,XS,x0,0);
figure();
subplot(2,1,1);
plot(T,XS(:,1));
ax=gca();
ax.XScale='log';

subplot(2,1,2);
plot(T,llh(:,1));
ax=gca();
ax.XScale='log';
hold('on');
plot(T,-.5*(log(2*pi)+log(V(:,1))),'DisplayName','LLH @Peak of gaussian');

[mle,mle_llh] = ou.estimateMLE(T,XS);