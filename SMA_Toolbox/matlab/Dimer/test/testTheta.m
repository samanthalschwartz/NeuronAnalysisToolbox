function testTheta(DA,DB,dt,r0,N)
    v=2*dt*(DA+DB);
    ds = repmat([r0;0],1,N) + randn(2,N)*sqrt(v);
    rs = sqrt(sum(ds.^2));
    thetas = atan2(ds(2,:),ds(1,:));
    
    f_r = @(r) (1/v)*exp( - (1/(2*v))* (r0^2 + r.^2)).*besseli(0, (r0/v)*r);
    figure();
    subplot(1,2,1);
    H=histogram(rs,'DisplayName','Sample','Normalization','probability');
    pdfs= .5*(f_r(H.BinEdges(2:end)) + f_r(H.BinEdges(1:end-1))) ./ diff(H.BinEdges);
    hold();
    plot(.5*(H.BinEdges(2:end) +H.BinEdges(1:end-1)), pdfs, 'k-');
    title('r distribution');
    xlabel('r');
    
end