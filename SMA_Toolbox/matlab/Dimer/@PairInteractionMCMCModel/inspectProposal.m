function inspectProposal(obj, Obs, P, PT, smc, obsIdx)
    %
    % Particle inspection
    %
    
    [Ps,weights,llh] = smc.getAllParticles(obsIdx);
    Ps=Ps{1};
    N = size(weights,1);
    llhAll = zeros(9,obsIdx,N);
    for n=1:N
        llhAll(:,:,n) = smc.computeLLH_debug(Obs(1:obsIdx,:)',Ps(:,:,n));
    end
    Nbound = sum(Ps(1,end,:));
    meColor = [0,0,0];
    p1Color = [1,0,0];
    p2Color = [0,0,1];
    bndColor = [0,1,0];
    tbndColor = [1,0,1];
    ms = 6;

    F1=figure('Position',[10,10,800,800]);
    cla();
    hold('on');
    
    for n=1:N
        Pn = Ps(:,:,n)';
        H=plot(Pn(obsIdx,2),Pn(obsIdx,3), 'Marker','v','MarkerEdgeColor', p1Color, 'MarkerFaceColor', p1Color*.5, 'MarkerSize',4);
        H.Annotation.LegendInformation.IconDisplayStyle='off';
        H=plot(Pn(obsIdx,4),Pn(obsIdx,5), 'Marker','^','MarkerEdgeColor', p2Color, 'MarkerFaceColor', p2Color*.5, 'MarkerSize',4);
        H.Annotation.LegendInformation.IconDisplayStyle='off';
        if Pn(obsIdx,1)
            vest = obj.transV(Pn(obsIdx,2:5));
            H=plot(Pn(obsIdx,[2,4]),Pn(obsIdx,[3,5]),'-','Color',[.5,1,.5],'LineWidth',0.5);
            H.Annotation.LegendInformation.IconDisplayStyle='off';
            H=plot(vest(1),vest(2),'+','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,.3],'MarkerSize',4);
            H.Annotation.LegendInformation.IconDisplayStyle='off';
       end
    
    end
%     if obsIdx>1
%             plot(Obs(obsIdx-1,2), Obs(obsIdx-1,3), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color*.5, 'MarkerSize',ms,'DisplayName','P1 last obs');
%             plot(Obs(obsIdx-1,4), Obs(obsIdx-1,5), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color*.5, 'MarkerSize',ms,'DisplayName','P2 last obs');
%             PairInteractionMCMCModel.plotSE(Obs(obsIdx-1,2:3), Obs(obsIdx-1,6:7), p1Color*0.5, 2, 'P1 last obs SE');
%             PairInteractionMCMCModel.plotSE(Obs(obsIdx-1,4:5), Obs(obsIdx-1,8:9), p2Color*0.5, 2, 'P2 last obs SE');
%     end
    plot(Obs(obsIdx,2), Obs(obsIdx,3), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 cur obs');
    plot(Obs(obsIdx,4), Obs(obsIdx,5), 'Marker','s','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 cur obs');
    PairInteractionMCMCModel.plotSE(Obs(obsIdx,2:3), Obs(obsIdx,6:7), p1Color, 2, 'P1 cur obs SE');
    PairInteractionMCMCModel.plotSE(Obs(obsIdx,4:5), Obs(obsIdx,8:9), p2Color, 2, 'P2 cur obs SE');

    if obsIdx>1
        plot(PT(obsIdx-1,2),PT(obsIdx-1,3), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color*.5, 'MarkerSize',ms,'DisplayName','[old] P1 true');
        plot(PT(obsIdx-1,4),PT(obsIdx-1,5), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color*.5, 'MarkerSize',ms,'DisplayName','[old] P2 true');
        plot(P(obsIdx-1,2),P(obsIdx-1,3), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color*.5, 'MarkerSize',ms,'DisplayName','[old] P1 est');
        plot(P(obsIdx-1,4),P(obsIdx-1,5), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color*.5, 'MarkerSize',ms,'DisplayName','[old] P2 est');
        if PT(obsIdx-1,1)
            vtrue = obj.transV(PT(obsIdx-1,2:5));
            plot(PT(obsIdx-1,[2,4]),PT(obsIdx-1,[3,5]),'-','Color',tbndColor*.5,'LineWidth',1,'DisplayName',sprintf('[old] True Bound R:%.3f Phi:%.3f',vtrue(3),vtrue(4)));
            plot(vtrue(1),vtrue(2),'d','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0]*.5,'DisplayName','[old] True center Position');
        end
        if P(obsIdx-1,1)
            vest = obj.transV(P(obsIdx-1,2:5));
            plot(P(obsIdx-1,[2,4]),P(obsIdx-1,[3,5]),'-','Color',bndColor,'LineWidth',1,'DisplayName',sprintf('[old] Est Bound R:%.3f Phi:%.3f',vest(3),vest(4)));
            plot(vest(1),vest(2),'o','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0]*.5,'DisplayName','[old] Est center Position');
        end
    end

    plot(PT(obsIdx,2),PT(obsIdx,3), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 true');
    plot(PT(obsIdx,4),PT(obsIdx,5), 'Marker','d','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 true');
    plot(P(obsIdx,2),P(obsIdx,3), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p1Color, 'MarkerSize',ms,'DisplayName','P1 est');
    plot(P(obsIdx,4),P(obsIdx,5), 'Marker','o','MarkerEdgeColor', meColor, 'MarkerFaceColor', p2Color, 'MarkerSize',ms,'DisplayName','P2 est');
    if PT(obsIdx,1)
        vtrue = obj.transV(PT(obsIdx,2:5));
        plot(PT(obsIdx,[2,4]),PT(obsIdx,[3,5]),'-','Color',tbndColor,'LineWidth',1.5,'DisplayName',sprintf('True Bound R:%.3f Phi:%.3f',vtrue(3),vtrue(4)));
        plot(vtrue(1),vtrue(2),'d','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0],'DisplayName','True center Position');
    end
    if P(obsIdx,1)
        vest = obj.transV(P(obsIdx,2:5));
        plot(P(obsIdx,[2,4]),P(obsIdx,[3,5]),'-','Color',bndColor,'LineWidth',1.5,'DisplayName',sprintf('Est Bound R:%.3f Phi:%.3f',vest(3),vest(4)));
        plot(vest(1),vest(2),'o','MarkerEdgeColor',meColor,'MarkerFaceColor',[1,1,0],'DisplayName','Est center Position');
    end
    

    legend('Location','BestOutside');
    axis('equal');
    title(sprintf('Observation %i/%i.  Bound %i/%i',obsIdx,size(Obs,1),Nbound,N ));
%         xlim([minX,maxX]);
%         ylim([minY,maxY]);

%     F2=figure('Position',[810,10,1100,800]);
    
end
