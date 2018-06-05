function test1DFree(y0,DA,DB,dt,N)
    y0=y0(:);
    vA=2*DA*dt;
    vB=2*DB*dt;

    gaussianLLH = @(d,v) -(size(d,2)/2*log(2*pi*v))-sum(d.^2,2)/(2*v);
    gaussianLLH2D = @(d,Sigma) -.5*(numel(d)*log(2*pi) + log(det(Sigma)) + d'*(Sigma\d));
    
    y1=[y0(1)+randn(1,N)*sqrt(vA); y0(2)+randn(1,N)*sqrt(vB)];
    pA1 = gaussianLLH(y1(1,:)'-y0(1),vA) ;
    pB1 =  gaussianLLH(y1(2,:)'-y0(2),vB);
    pY1 = pA1+pB1;
    
    A = [ .5 .5; 1 -1];
    z0 = A*y0;
    z1 = A*y1;
    vZ = A* [vA, 0;0, vB]*A';
    pZ1 = arrayfun(@(za,zb) gaussianLLH2D([za;zb]-z0,vZ),z1(1,:), z1(2,:));
    

    figure();
    subplot(1,2,1);
    plot(y0(1),y0(2),'ko','DisplayName','y0');
    hold();
    plot(y1(1,:),y1(2,:),'r.','DisplayName','y1');
    axis('equal');
    legend('location','best');
    xlabel('Particle A');
    ylabel('Particle B');

    subplot(1,2,2);
    plot(1:N,pA1,'r.','DisplayName','p(A1)');
    hold();
    plot(1:N,pB1,'b.','DisplayName','p(B1)');
    plot(1:N,pY1,'k.','DisplayName','p(Y1)=p(A1)p(B1)');
%     plot(1:N,pC1,'ro','DisplayName','p(C1)');
%     plot(1:N,pD1,'bo','DisplayName','p(D1)');
    plot(1:N,pZ1,'ko','DisplayName','p(Z1)=p(C1)p(D1)');
    legend('location','best');
    xlabel('Sequence');
    ylabel('LLH');
end
