function plotNormalPDF(N)

    sigma=logspace(-6,6,100);
    
    expected_LLH = arrayfun(@(s) mean(log(normpdf(s*randn(N,1),0,s))), sigma);
    
    f=figure();
    semilogx(sigma, expected_LLH,'-r');
    xlabel('sigma');
    ylabel('$\langle \mathrm{LLH} \rangle$','interpreter','latex');


end