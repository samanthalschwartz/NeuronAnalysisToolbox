function testJumpLikelihoods(DA,DB,DAB,Dphi,sigma_r,rho,dt,N)
    %
    % Spaces:
    %  Y-space: [a_x a_y b_x b_y]
    %  U-space: [c_x c_y d_x d_y]
    %  V-space: [c_x c_y r phi]
    %

    % U Transform
    Umat=[.5 0 .5 0; 0 .5 0 .5; 1 0 -1 0; 0 1 0 -1]; % Linear transform to U space
    Umatinv = inv(Umat);
    U = @(ys) Umat*ys;
    Uinv = @(ys) Umatinv*ys;

    % V Transform
    V = @(ys) [.5*sum(ys([1,3],:)); .5*sum(ys([2,4],:)); sqrt(sum((ys(1:2,:)-ys(3:4,:)).^2)); atan2(ys(2,:)-ys(4,:),ys(1,:)-ys(3,:))];
    Vinv = @(vs) [vs(1,:)+.5*vs(3,:).*cos(vs(4,:)); vs(2,:)+.5*vs(3,:).*sin(vs(4,:)); vs(1,:)-.5*vs(3,:).*cos(vs(4,:)); vs(2,:)-.5*vs(3,:).*sin(vs(4,:))];    

    %Initial conditions
    y0=[rho,rho,0,rho]'; % Initial positions in Y space
    u0=U(y0); % Initial positions in U space
    v0=V(y0); % Initial positions in V space

    % Definition of provability distribution for single step
    SigmaFree_Y = 2*dt*diag([DA DA DB DB]); % Variance for free case in Y space
    SigmaBound_V = diag([2*dt*DAB 2*dt*DAB sigma_r^2 2*dt*Dphi]); % Variance for bound case in V space
    
    %Simulate a single step
    y1_F = repmat(y0,1,N) + chol(SigmaFree_Y)*randn(4,N); % chol is effectively square root
    v1_B = repmat(v0,1,N) + chol(SigmaBound_V)*randn(4,N); % chol is effectively square root

    %Transform simulations to other coordinates
    u1_F = U(y1_F);
    v1_F = V(y1_F);

    y1_B = Vinv(v1_B);
    u1_B = U(y1_B);
   

    %Compute likelihoods
    SigmaFree_U = Umat*SigmaFree_Y*Umat';

    LLH_F.y1_F = sum(multiGaussianLLH(repmat(y0,1,N) - y1_F, SigmaFree_Y));
    LLH_F.u1_F = sum(multiGaussianLLH(repmat(u0,1,N) - u1_F, SigmaFree_U));
    LLH_F.v1_F = sum(multiGaussianLLH(repmat(Vinv(v0),1,N) - Vinv(v1_F), SigmaFree_Y)) + log(v1_F(3,:)); 
    LLH_F.y1_B = sum(multiGaussianLLH(repmat(y0,1,N) - y1_B, SigmaFree_Y));
    LLH_F.u1_B = sum(multiGaussianLLH(repmat(u0,1,N) - u1_B, SigmaFree_U));
    LLH_F.v1_B = sum(multiGaussianLLH(repmat(Vinv(v0),1,N) - Vinv(v1_B), SigmaFree_Y)) + log(v1_B(3,:)); 

    LLH_B.v1_F = sum(multiGaussianLLH(repmat(v0,1,N) - v1_F, SigmaBound_V));
    LLH_B.v1_B = sum(multiGaussianLLH(repmat(v0,1,N) - v1_B, SigmaBound_V));
    LLH_B.y1_F = sum(multiGaussianLLH(V(repmat(y0,1,N)) - V(y1_F), SigmaBound_V)) - log(v1_F(3,:));
    LLH_B.y1_B = sum(multiGaussianLLH(V(repmat(y0,1,N)) - V(y1_B), SigmaBound_V)) - log(v1_B(3,:));
    LLH_B.u1_F = sum(multiGaussianLLH(V(Uinv(repmat(u0,1,N))) - V(Uinv(u1_F)), SigmaBound_V)) - log(v1_F(3,:));
    LLH_B.u1_B = sum(multiGaussianLLH(V(Uinv(repmat(u0,1,N))) - V(Uinv(u1_B)), SigmaBound_V)) - log(v1_B(3,:));

%     Marg_PDF_F.phi = @(phis) multiGaussianLLH(
% 
%     Marg_PDF_B.phi = @(phis) exp(multiGaussianLLH(repmat(v0(4),1,numel(phis)) - phis, SigmaBound_V(4,4)));

    %Plot positions
    ms=6;
    figure('Position',[10,10,1200,1000]);
    ax_pos_F=subplot(2,2,1);
    plot(y0(1),y0(2),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[.5,0,0],'DisplayName','$a_0$');
    hold();
    plot(y0(3),y0(4),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,.5],'DisplayName','$b_0$');
    plot(u0(1),u0(2),'o','MarkerFaceColor',[1,.5,0],'MarkerEdgeColor',[.5,.5,0],'DisplayName','$c_0 =(a_0+b_0)/2$');
    plot(u0(3),u0(4),'o','MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,1,.5],'DisplayName','$d_0= a_0-b_0$');
    plot(u0(1),u0(2),'o','MarkerFaceColor',[1,.5,0],'MarkerEdgeColor',[.5,.5,0],'DisplayName','$c_0 =(a_0+b_0)/2$');
    plot(u0(3),u0(4),'o','MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,1,.5],'DisplayName','$d_0= a_0-b_0$');
    plot(y1_F(1,:), y1_F(2,:),'.','MarkerSize',ms,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'DisplayName','Free: $a_1$');
    plot(y1_F(3,:), y1_F(4,:),'.','MarkerSize',ms,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'DisplayName','Free: $b_1$');
    plot(u1_F(1,:), u1_F(2,:),'.','MarkerSize',ms,'MarkerFaceColor',[1,.5,0],'MarkerEdgeColor',[1,.5,0],'DisplayName','Free: $c_1=(a_1+b_1)/2$');
    plot(u1_F(3,:), u1_F(4,:),'.','MarkerSize',ms,'MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,1,1],'DisplayName','Free: $d_1=a_1-b_1$');

    H=legend('Location','best');
    H.Interpreter='latex';
    axis('equal');
    title('Free');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');

    ax_histr_F=subplot(4,2,2);
    histogram(v1_F(3,:));
    xlabel('Distance: r (um)','interpreter','latex');
    title('Free','interpreter','latex');

    ax_histphi_F=subplot(4,2,4);
    histogram(v1_F(4,:));
    xlabel('Angle: $\phi$ (rad)','interpreter','latex');
    title('Free','interpreter','latex');


    ax_pos_B=subplot(2,2,3);
    plot(y0(1),y0(2),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[.5,0,0],'DisplayName','$a_0$');
    hold();
    plot(y0(3),y0(4),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,.5],'DisplayName','$b_0$');
    plot(u0(1),u0(2),'o','MarkerFaceColor',[1,.5,0],'MarkerEdgeColor',[.5,.5,0],'DisplayName','$c_0 =(a_0+b_0)/2$');
    plot(u0(3),u0(4),'o','MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,1,.5],'DisplayName','$d_0= a_0-b_0$');

    plot(y1_B(1,:), y1_B(2,:),'.','MarkerSize',ms,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'DisplayName','Bound: $a_1$');
    plot(y1_B(3,:), y1_B(4,:),'.','MarkerSize',ms,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'DisplayName','Bound: $b_1$');
    plot(u1_B(1,:), u1_B(2,:),'.','MarkerSize',ms,'MarkerFaceColor',[1,.5,0],'MarkerEdgeColor',[1,.5,0],'DisplayName','Bound: $c_1=(a_1+b_1)/2$');
    plot(u1_B(3,:), u1_B(4,:),'.','MarkerSize',ms,'MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,1,1],'DisplayName','Bound: $d_1=a_1-b_1$');
%     plot(v1_B(3,:).*cos(v1_B(4,:)), v1_B(3,:).*sin(v1_B(4,:)),'o','MarkerSize',ms,'MarkerFaceColor','none','MarkerEdgeColor',[0,1,1],'DisplayName','Bound: $d_1=a_1-b_1$');

    
    H=legend('Location','best');
    H.Interpreter='latex';
    axis('equal');
    title('Bound');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');

    ax_histr_B=subplot(4,2,6);
    H=histogram(v1_B(3,:),'Normalization','probability','DisplayName','Simulated sample','BinMethod','scott');
    Ds=[H.BinEdges(1:end-1);H.BinEdges(2:end)];
    pdfs = gaussianPDF_interval(Ds-rho,SigmaBound_V(3,3));
    hold();
    plot(mean(Ds),pdfs,'k-','DisplayName','Marginal Pdf (per bin)');
    legend('location','best');
    xlabel('Distance: r (um)','interpreter','latex');
    title('Bound','interpreter','latex');

    ax_histphi_B=subplot(4,2,8);
    H=histogram(v1_B(4,:),'Normalization','probability','DisplayName','Simulated sample','BinMethod','scott');
    Ds=[H.BinEdges(1:end-1);H.BinEdges(2:end)];
    pdfs = gaussianPDF_interval(Ds-v0(4),SigmaBound_V(4,4));
    hold();
    plot(mean(Ds),pdfs,'k-','DisplayName','Marginal Pdf (per bin)');
    legend('location','best');
    xlabel('Angle: $\phi$ (rad)','interpreter','latex');
    title('Bound','interpreter','latex');

    linkaxes([ax_pos_F,ax_pos_B],'xy');
%     linkaxes([ax_histr_F, ax_histr_B],'xy');
%     linkaxes([ax_histphi_F,ax_histphi_B],'xy');

    %Plot LLHs
    
    figure();
    subplot(2,2,1);
    plot(1:N,LLH_F.y1_F,'-','Color',[1,.5,0],'DisplayName','llh(y1_F|y0, Z_1=F)');
    hold();
    plot(1:N,LLH_F.u1_F,'.','Color',[1,.5,0],'DisplayName','llh(u1_F|u0, Z_1=F)');
    plot(1:N,LLH_B.y1_F,'-','Color',[0,.5,1],'DisplayName','llh(y1_F|v0, Z_1=B)');
    plot(1:N,LLH_B.u1_F,'.','Color',[0,.5,1],'DisplayName','llh(u1_F|u0, Z_1=B)');
    title('TrueState:Free [Y/U-space]');
    xlabel('Sample');
    ylabel('Log Likelihood');
    legend('Location','best');

    subplot(2,2,2);
    plot(1:N,LLH_F.v1_F,'-','Color',[1,.5,0],'DisplayName','llh(v1_F|v0, Z_1=F)');
    hold();
    plot(1:N,LLH_B.v1_F,'-','Color',[0,.5,1],'DisplayName','llh(v1_F|v0, Z_1=B)');
    title('TrueState:Free [V-space]');
    xlabel('Sample');
    ylabel('Log Likelihood');
    legend('Location','best');

    subplot(2,2,3);
    plot(1:N,LLH_F.y1_B,'-','Color',[.5,.25,0],'DisplayName','llh(y1_B|y0, Z_1=F)');
    hold();
    plot(1:N,LLH_F.u1_B,'.','Color',[.5,.25,0],'DisplayName','llh(u1_B|u0, Z_1=F)');
    plot(1:N,LLH_B.y1_B,'-','Color',[0,.25,.5],'DisplayName','llh(y1_B|v0, Z_1=B)');
    plot(1:N,LLH_B.u1_B,'.','Color',[0,.25,.5],'DisplayName','llh(u1_B|u0, Z_1=B)');
    title('TrueState:Bound [Y/U-space]');
    xlabel('Sample');
    ylabel('Log Likelihood');
    legend('Location','best');

    subplot(2,2,4);
    plot(1:N,LLH_F.v1_B,'-','Color',[.5,.25,0],'DisplayName','llh(v1_B|v0, Z_1=F)');
    hold();
    plot(1:N,LLH_B.v1_B,'-','Color',[0,.25,.5],'DisplayName','llh(v1_B|v0, Z_1=B)');
    title('TrueState:Bound [V-space]');
    xlabel('Sample');
    ylabel('Log Likelihood');
    legend('Location','best');



end

function cdfs= gaussianCDF(ds,V)
    % ds - size[1,N] distances
    % V - variance
    cdfs = .5*(1+erf(ds./sqrt(2*V)));
end

function pdfs=gaussianPDF_interval(ds,V)
    % ds - size[2,N] ds(1,:) interval start, ds(2,:) interval end Note: all(ds(1,:) < ds(2,:))
    % V - variance scalar
    pdfs = gaussianCDF(ds(2,:),V)-gaussianCDF(ds(1,:),V);
end


function llhs=multiGaussianLLH(ds,Sigma)
    % [in]
    % ds: size:[k,N] each row is a distance vector in R^k
    % SigmaL size:[k,k]
    % [out]
    % llhs: size:[1,N] Log-likelihoods
    llhs=-.5*(size(ds,1)*log(2*pi)+log(det(Sigma))) - .5* ds.*(Sigma\ds);
end

