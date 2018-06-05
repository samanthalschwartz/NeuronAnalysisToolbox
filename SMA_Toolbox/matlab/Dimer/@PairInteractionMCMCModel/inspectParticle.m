function inspectParticle(obj, data_idx, P, smc)
    %
    % Particle inspection
    %

    F1=figure('Position',[10,10,800,800]);
    obsAx = gca();
    F2=figure('Position',[810,10,1100,800]);
    D=obj.Data(data_idx);
    Obs = D.obs;
    obsIdx = 1;
    PT = D.particleTrue;
    maxX = max([max(Obs(:,[2,4])+3*Obs(:,[6,8])), max(P(:,[2,4])), max(PT(:,[2,4]))]);
    minX = min([min(Obs(:,[2,4])-3*Obs(:,[6,8])), min(P(:,[2,4])), min(PT(:,[2,4]))]);
    maxY = max([max(Obs(:,[3,5])+3*Obs(:,[7,9])), max(P(:,[3,5])), max(PT(:,[3,5]))]);
    minY = min([min(Obs(:,[3,5])-3*Obs(:,[7,9])), min(P(:,[3,5])), min(PT(:,[3,5]))]);

    Pllh = smc.computeLLH_debug(D.obs',P')';
    Pllh_total = sum(Pllh,2);
    PTllh = smc.computeLLH_debug(D.obs',PT')';
    PTllh_total = sum(PTllh,2);
    Nobs = D.Nobs;
    controlH = 0.035; % Contols height
    uicontrol('Parent',F1,'Style','slider','Min',1, 'Max',Nobs,'Value',1,'SliderStep',[1/(Nobs-1),1/(Nobs-1)],...
                  'Units','normalized','Position',[0.0, 0, 0.7, controlH], 'Callback',@slider_CB);
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.7, 0, 0.1, controlH],'String','[I] Proposal', 'Callback', @inspectProposal_CB);
    

    meColor = [0,0,0];
    p1Color = [1,0,0];
    p2Color = [0,0,1];
    bndColor = [0,1,0];
    tbndColor = [1,0,1];
    ms = 6;

    plotObs();

    function slider_CB(hObj,~)
        k=round(hObj.Value);
        obsIdx=k;
        plotObs();
    end

    function inspectProposal_CB(~,~)
        obj.inspectProposal(Obs,P,PT,smc,obsIdx);
    end
    
    function plotObs()
        axes(obsAx);
        cla();
        hold('on');
        if obsIdx>1
            plot(Obs(obsIdx-1,2), Obs(obsIdx-1,3), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color*.5, 'MarkerSize',ms,'DisplayName','P1 last obs');
            plot(Obs(obsIdx-1,4), Obs(obsIdx-1,5), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color*.5, 'MarkerSize',ms,'DisplayName','P2 last obs');
            PairInteractionMCMCModel.plotSE(Obs(obsIdx-1,2:3), Obs(obsIdx-1,6:7), p1Color*0.5, 2, 'P1 last obs SE');
            PairInteractionMCMCModel.plotSE(Obs(obsIdx-1,4:5), Obs(obsIdx-1,8:9), p2Color*0.5, 2, 'P2 last obs SE');
        end
        plot(Obs(obsIdx,2), Obs(obsIdx,3), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 cur obs');
        plot(Obs(obsIdx,4), Obs(obsIdx,5), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 cur obs');
        PairInteractionMCMCModel.plotSE(Obs(obsIdx,2:3), Obs(obsIdx,6:7), p1Color, 2, 'P1 cur obs SE');
        PairInteractionMCMCModel.plotSE(Obs(obsIdx,4:5), Obs(obsIdx,8:9), p2Color, 2, 'P2 cur obs SE');
        plot(PT(obsIdx,2),PT(obsIdx,3), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 true');
        plot(PT(obsIdx,4),PT(obsIdx,5), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 true');

        plot(P(obsIdx,2),P(obsIdx,3), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 est');
        plot(P(obsIdx,4),P(obsIdx,5), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 est');
        if PT(obsIdx,1)
            vtrue = obj.transV(PT(obsIdx,2:5));
            plot(PT(obsIdx,[2,4]),PT(obsIdx,[3,5]),'-','Color',tbndColor,'LineWidth',1.5,'DisplayName',sprintf('True Bound R:%.3f Phi:%.3f',vtrue(3),vtrue(4)));
            plot(vtrue(1),vtrue(2),'o','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0]*0.5,'DisplayName','True center Position');
        end
        if P(obsIdx,1)
            vest = obj.transV(P(obsIdx,2:5));
            plot(P(obsIdx,[2,4]),P(obsIdx,[3,5]),'-','Color',bndColor,'LineWidth',1.5,'DisplayName',sprintf('Est Bound R:%.3f Phi:%.3f',vest(3),vest(4)));
            plot(vest(1),vest(2),'o','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0],'DisplayName','Est center Position');
        end
        legend('Location','BestOutside');
        xlim([minX,maxX]);
        ylim([minY,maxY]);

        axis('equal');
        title(sprintf('Observation %i/%i',obsIdx,Nobs));

        figure(F2);
        subplot(3,5,1);
        cla();
        hold('on');
        stairs(P(:,1),'r-','DisplayName','Particle State');
        stairs(PT(:,1),'k-','DisplayName','True State');
        plot(obsIdx,P(obsIdx,1),'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle State');
        plot(obsIdx,PT(obsIdx,1),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True State');
        legend('Location','Best');
        xlim([0,Nobs]);
        title('State');

        subplot(3,5,2);
        cla();
        hold('on');
        plot(Pllh_total,'r-','DisplayName','Particle LLH');
        plot(PTllh_total,'k-','DisplayName','True LLH');
        plot(obsIdx,Pllh_total(obsIdx),'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'DisplayName',sprintf('Particle LLH: %.5g',Pllh_total(data_idx)));
        plot(obsIdx,PTllh_total(obsIdx),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName',sprintf('True LLH: %.5g',PTllh_total(data_idx)));
        legend('Location','Best');
        xlim([0,Nobs]);
        yl=[min([Pllh_total(:); PTllh_total(:); Pllh(:); PTllh(:)]),max([Pllh_total(:); PTllh_total(:); Pllh(:); PTllh(:)])];
        ylim(yl);
        title('Overall LLH');

        subplot(3,5,3);
        cla();
        hold('on');
        plot(Pllh(:,1),'r-','DisplayName','Particle Z LLH');
        plot(PTllh(:,1),'k-','DisplayName','True Z LLH');
        plot(obsIdx,Pllh(obsIdx,1),'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle Z LLH');
        plot(obsIdx,PTllh(obsIdx,1),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True Z LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Z LLH');

        subplot(3,5,6);
        cla();
        hold('on');
        plot(Pllh(:,2),'b-','DisplayName','Particle y1 LLH');
        plot(PTllh(:,2),'k-','DisplayName','True y1 LLH');
        plot(obsIdx,Pllh(obsIdx,2),'o','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',5,'DisplayName','Particle y1 LLH');
        plot(obsIdx,PTllh(obsIdx,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True y1 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Y1 LLH');

        subplot(3,5,7);
        cla();
        hold('on');
        plot(Pllh(:,3),'b-','DisplayName','Particle y2 LLH');
        plot(PTllh(:,3),'k-','DisplayName','True y2 LLH');
        plot(obsIdx,Pllh(obsIdx,3),'o','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',5,'DisplayName','Particle y2 LLH');
        plot(obsIdx,PTllh(obsIdx,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True y2 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Y2 LLH');

        subplot(3,5,8);
        cla();
        hold('on');
        plot(Pllh(:,4),'b-','DisplayName','Particle y3 LLH');
        plot(PTllh(:,4),'k-','DisplayName','True y3 LLH');
        plot(obsIdx,Pllh(obsIdx,4),'o','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',5,'DisplayName','Particle y3 LLH');
        plot(obsIdx,PTllh(obsIdx,4),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True y3 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Y3 LLH');

        subplot(3,5,9);
        cla();
        hold('on');
        plot(Pllh(:,5),'b-','DisplayName','Particle y4 LLH');
        plot(PTllh(:,5),'k-','DisplayName','True y4 LLH');
        plot(obsIdx,Pllh(obsIdx,5),'o','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',5,'DisplayName','Particle y4 LLH');
        plot(obsIdx,PTllh(obsIdx,5),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True y4 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Y4 LLH');

        subplot(3,5,10);
        cla();
        hold('on');
        plot(sum(Pllh(:,2:5),2),'b-','DisplayName','Particle Y LLH');
        plot(sum(PTllh(:,2:5),2),'k-','DisplayName','True Y LLH');
        plot(obsIdx,sum(Pllh(obsIdx,2:5),2),'o','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',5,'DisplayName','Particle Y LLH');
        plot(obsIdx,sum(PTllh(obsIdx,2:5),2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True Y LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Ytotal LLH');
        
        subplot(3,5,11);
        cla();
        hold('on');
        plot(Pllh(:,6),'m-','DisplayName','Particle x1 LLH');
        plot(PTllh(:,6),'k-','DisplayName','True x1 LLH');
        plot(obsIdx,Pllh(obsIdx,6),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle x1 LLH');
        plot(obsIdx,PTllh(obsIdx,6),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True x1 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('X1 LLH');

        subplot(3,5,12);
        cla();
        hold('on');
        plot(Pllh(:,7),'m-','DisplayName','Particle x2 LLH');
        plot(PTllh(:,7),'k-','DisplayName','True x2 LLH');
        plot(obsIdx,Pllh(obsIdx,7),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle x2 LLH');
        plot(obsIdx,PTllh(obsIdx,7),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True x2 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('X2 LLH');

        subplot(3,5,13);
        cla();
        hold('on');
        plot(Pllh(:,8),'m-','DisplayName','Particle x3 LLH');
        plot(PTllh(:,8),'k-','DisplayName','True x3 LLH');
        plot(obsIdx,Pllh(obsIdx,8),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle x3 LLH');
        plot(obsIdx,PTllh(obsIdx,8),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True x3 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('X3 LLH');

        subplot(3,5,14);
        cla();
        hold('on');
        plot(Pllh(:,9),'m-','DisplayName','Particle x4 LLH');
        plot(PTllh(:,9),'k-','DisplayName','True x4 LLH');
        plot(obsIdx,Pllh(obsIdx,9),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle x4 LLH');
        plot(obsIdx,PTllh(obsIdx,9),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True x4 LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('X4 LLH');

        subplot(3,5,15);
        cla();
        hold('on');
        plot(sum(Pllh(:,6:9),2),'m-','DisplayName','Particle X LLH');
        plot(sum(PTllh(:,6:9),2),'k-','DisplayName','True X LLH');
        plot(obsIdx,sum(Pllh(obsIdx,6:9),2),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','Particle X LLH');
        plot(obsIdx,sum(PTllh(obsIdx,6:9),2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'DisplayName','True X LLH');
        legend('Location','Best');
        xlim([0,Nobs]);
        ylim(yl);
        title('Xtotal LLH');
    end  
end
