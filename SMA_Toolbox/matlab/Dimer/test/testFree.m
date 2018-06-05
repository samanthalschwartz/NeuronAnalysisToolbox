function testFree(a0,b0,DA,DB,dt,N)
    c0=(a0+b0)/2;
    d0=a0-b0;
    sigmaA=sqrt(2*DA*dt);
    sigmaB=sqrt(2*DB*dt);
    a1=repmat(a0,N,1)+randn(N,2).*sigmaA;
    b1=repmat(b0,N,1)+randn(N,2).*sigmaB;
    c1 = (a1+b1)/2;
    d1 = a1-b1;
    
    var_a1 = 2*DA*dt;
    var_b1 = 2*DB*dt;
    var_c1 = .25*(var_a1+var_b1);
    var_d1 = 2*(DA+DB)*dt;
    
    gaussianLLH = @(d,v) -(size(d,2)/2*log(2*pi*v))-sum(d.^2,2)/(2*v);
    
    LLH_F_a = gaussianLLH(a1,var_a1);    
    LLH_F_b = gaussianLLH(b1,var_b1);
    LLH_F_c = gaussianLLH(c1,var_c1);
    LLH_F_d = gaussianLLH(d1,var_d1);
    LLH_F_ab = LLH_F_a+LLH_F_b;
    LLH_F_cd = LLH_F_c+LLH_F_d;
    
    
    figure('Position',[10,10,1400,600]);
    subplot(1,2,1);
    plot(a0(1),a0(2),'ro','DisplayName','a_0','MarkerFaceColor',[.5,0,0]);
    hold();
    plot(b0(1),b0(2),'bo','DisplayName','b_0','MarkerFaceColor',[0,0,.5]);
    plot(a1(:,1),a1(:,2),'r.','DisplayName','a_1');
    plot(b1(:,1),b1(:,2),'b.','DisplayName','b_1');
    plot(c0(:,1),c0(:,2),'go','DisplayName','c_0','MarkerFaceColor',[0,.5,0]);
    plot(c1(:,1),c1(:,2),'g.','DisplayName','c_1');
    plot(d0(:,1),d0(:,2),'mo','DisplayName','d_0','MarkerFaceColor',[.5,0,.5]);
    plot(d1(:,1),d1(:,2),'m.','DisplayName','d_1');
    legend('location','best');
    axis('equal');
    subplot(1,2,2);
    plot(1:N, LLH_F_a,'r-','DisplayName','LLH(a_i)');
    hold();
    plot(1:N, LLH_F_b,'b-','DisplayName','LLH(b_i)');
    plot(1:N, LLH_F_c,'g-','DisplayName','LLH(c_i)');
    plot(1:N, LLH_F_d,'m-','DisplayName','LLH(d_i)');
    plot(1:N, LLH_F_cd,'k-','LineWidth',2,'DisplayName','LLH(c_i, d_i)');
    plot(1:N, LLH_F_ab,'y-','LineWidth',2,'DisplayName','LLH(a_i, b_i)');
    legend('location','best');
    xlabel('Sample');
    ylabel('LLH');
end