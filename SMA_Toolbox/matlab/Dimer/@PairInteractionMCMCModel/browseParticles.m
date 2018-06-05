function browseParticles(obj, exper_idx, Nparticles, ProposalType)
    %
    % Particle browser
    %
    if nargin<2
        exper_idx=1;
    end
    if nargin<3
        Nparticles=250;
    end
    if nargin<4
        ProposalType='Transition';
    end
    F1=figure('Position',[10,10,800,700]);
    F2=figure('Position',[810,10,800,700]);
    whitebg(F1);
%             whitebg(F2);
    F1.Color=[0 0 0];
%             F2.Color=[0 0 0];


    [smc,D] = obj.makeExperimentSMC(exper_idx);
    Ps = [];
    PsLLH = [];
    Os = {D(:).obs}';
    PTs = {D(:).particleTrue}';
    Ndata = numel(Os);




    % Find the min and max positions and times over all trajectories: observed trajs and particle trajs
    allTs = [Os;PTs];
    xmin = min(cellfun(@(T) min(min(T(:,[2,4]))), allTs));
    xmax = max(cellfun(@(T) max(max(T(:,[2,4]))), allTs));
    ymin = min(cellfun(@(T) min(min(T(:,[3,5]))), allTs));
    ymax = max(cellfun(@(T) max(max(T(:,[3,5]))), allTs));
    tmin = min(cellfun(@(T) min(T(:,1)), Os));
    tmax = max(cellfun(@(T) max(T(:,1)), Os));
    figure(F1);
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    zlim([tmin,tmax]);
    view(-45,30);
    asp(1) = (xmax-xmin)/(ymax-ymin);
    asp(2) = 1;
    asp(3) = mean([asp(1),asp(2)]);
    pbaspect(asp);


    bndAlpha=0.75;
    H=struct(); % Handles
    drawBonds=true; %Change this to control bond drawing

    figure(F1);
    grid('on');
    grid('minor');
    box('on');

    data_idx=1;

    resample();
    plotT();

    %Draw controls
    controlH = 0.035; % Contols height
    if Ndata>1
        uicontrol('Parent',F1,'Style','slider','Min',1, 'Max',Ndata,'Value',1,'SliderStep',[1/(Ndata-1),1/(Ndata-1)],...
                  'Units','normalized','Position',[0.0, 0, 0.5, controlH], 'Callback',@slider_CB);
    end
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.5, 0, 0.1, controlH],'String','Inspect', 'Callback',@inspect_CB);
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.6, 0, 0.1, controlH],'String','Resamp. T', 'Callback',@resampleT_CB);
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.7, 0, 0.1, controlH],'String','Resamp. O', 'Callback',@resampleO_CB);
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.8, 0, 0.1, controlH],'String','Resamp. C', 'Callback',@resampleC_CB);
    uicontrol('Parent',F1,'Style','pushbutton','Units','normalized','Position',[.9, 0, 0.1, controlH],'String','Resamp. M', 'Callback',@resampleM_CB);
    
    figure(F1);
    title(sprintf('Trajectory %i',1));
    xlabel('X (um)');
    ylabel('Y (um)');
    zlabel('T (s)');
    ax=gca();
    ax.TickDir='out';
    ax.BoxStyle='Full';
    legend('location','best');

    function slider_CB(hObj,~)
        k=round(hObj.Value);
        data_idx=k;
        plotT();
        title(sprintf('Trajectory %i',k));
    end

    function inspect_CB(varargin) 
        obj.inspectParticle(data_idx, Ps{data_idx}, smc);
    end

    function resampleT_CB(varargin)
        ProposalType='Transition';
        resample()
        plotT();
    end
    function resampleO_CB(varargin)
        ProposalType='Observation';
        resample()
        plotT();
    end
    function resampleC_CB(varargin)
        ProposalType='Combined';
        resample()
        plotT();
    end

    function resampleM_CB(varargin)
        exper_data_idx = obj.getExperimentDataIdxs(exper_idx);
        [Ps,weights,PsLLH] = obj.runParticleFilter(exper_data_idx, Nparticles, false);
        
        idx = randsample(Nparticles,1,true,weights);
        Ps = {Ps(:,:,idx)};
        PsLLH = PsLLH(idx);
        plotT();
    end

    function resample()
        smc.runParticleFilter(Nparticles, ProposalType);
        [Ps,PsLLH]=smc.sampleParticle();
        Ps = cellmap(@transpose,Ps);
    end

    function clearHandles(Hs)
        fs=fieldnames(Hs);
        N=numel(fs);
        for n=1:N
            F = Hs.(fs{n});
            if ishandle(F)
                delete(F)
            elseif iscell(F)
                for k=1:numel(F)
                    K = F{k};
                    if ishandle(K)
                        delete(K)
                    end
                end
            end
        end
    end

    function plotT()
        clearHandles(H);
        O = Os{data_idx};
        P = Ps{data_idx};
        PT = PTs{data_idx};
        N =size(P,1);
        times = O(:,[1,1]);
        figure(F1);
        hold('on');
        H.a=surface(O(:,[2,2]),O(:,[3,3]),times,'FaceColor',[0,0,0],'EdgeColor',[1,0,0],'LineWidth',1,'DisplayName','Obs A');
        H.b=surface(O(:,[4,4]),O(:,[5,5]),times,'FaceColor',[0,0,0],'EdgeColor',[0,0,1],'LineWidth',1,'DisplayName','Obs B');
        H.pta=surface(PT(:,[2,2]),PT(:,[3,3]),times,'FaceColor',[0,0,0],'EdgeColor',[1,0,0],'LineWidth',.5,'LineStyle','--','DisplayName','True A');
        H.ptb=surface(PT(:,[4,4]),PT(:,[5,5]),times,'FaceColor',[0,0,0],'EdgeColor',[0,0,1],'LineWidth',.5,'LineStyle','--','DisplayName','True B');
        bndTs=find(PT(:,1));
        if drawBonds            
            H.bds=cell(numel(bndTs),1);
            for i=1:numel(bndTs)
                b = PT(bndTs(i),:);
                t = O(bndTs(i),1);
                H.bds{i}=surface(b([2 2;4 4]),b([3 3;5 5]),t([1 1;1 1]),...
                                 'EdgeColor',[0,1,0],'EdgeAlpha',bndAlpha,'LineStyle','--');
                H.bds{i}.Annotation.LegendInformation.IconDisplayStyle='off';
            end
        end

        H.pa=surface(P(:,[2,2]),P(:,[3,3]),times,'FaceColor',[0,0,0],'EdgeColor',[1,.5,0],'DisplayName','Sampled A');
        H.pb=surface(P(:,[4,4]),P(:,[5,5]),times,'FaceColor',[0,0,0],'EdgeColor',[0,.5,1],'DisplayName','Sampled B');
        if drawBonds            
            PbndTs=find(P(:,1));
            H.Pbds=cell(N,1);
            for i=1:N
                if P(i,1) && ~PT(i,1)
                    C=[1,0,1];
                    LS='--';
                elseif ~P(i,1) && PT(i,1)
                    C=[1,0,1];
                    LS='-';
                elseif P(i,1) && PT(i,1)
                    C=[0,1,0];
                    LS='-';
                else
                    continue;
                end
                b = P(i,:);
                t = O(i,1);
                H.Pbds{i}=surface(b([2 2;4 4]),b([3 3;5 5]),t([1 1;1 1]),'LineStyle',LS,...
                                 'EdgeColor',C,'EdgeAlpha',bndAlpha); 
                H.Pbds{i}.Annotation.LegendInformation.IconDisplayStyle='off';
            end
        end

        figure(F2);
        subplot(2,2,1);
        cla();
        deltaObs = abs(O(:,2:5)-P(:,2:5));
        deltaTrue = abs(PT(:,2:5)-P(:,2:5));
        hold('on');
        plot(1:N,deltaTrue(:,1),'-','Color',[1,0,0],'DisplayName','\delta a_x');
        plot(1:N,deltaTrue(:,2),'--','Color',[1,0,0],'DisplayName','\delta a_y');
        plot(1:N,deltaTrue(:,3),'-','Color',[0,0,1],'DisplayName','\delta b_x');
        plot(1:N,deltaTrue(:,4),'--','Color',[0,0,1],'DisplayName','\delta b_y');
        plot(1:N,deltaObs(:,1),'-','Color',[1,0,1],'DisplayName','\delta a_x (obs)');
        plot(1:N,deltaObs(:,2),'--','Color',[1,0,1],'DisplayName','\delta a_y (obs)');
        plot(1:N,deltaObs(:,3),'-','Color',[0,1,1],'DisplayName','\delta b_x (obs)');
        plot(1:N,deltaObs(:,4),'--','Color',[0,1,1],'DisplayName','\delta b_y (obs)');
        plot(1:N,O(:,6),'-','Color',[0,0,0],'DisplayName','w a_x');
        plot(1:N,O(:,7),'--','Color',[0,0,0],'DisplayName','w a_y');
        plot(1:N,O(:,8),'-','Color',[.2,.2,.2],'DisplayName','w b_x');
        plot(1:N,O(:,9),'--','Color',[.2,.2,.2],'DisplayName','w b_y');
        title('Postion Estimation Accuracy');
        legend('location','best');
        ylabel('Absolute distance');
        subplot(2,2,2);
        cla();
        hold('on');
        distP = sqrt(sum((P(:,2:3)-P(:,4:5)).^2,2));
        distPT = sqrt(sum((PT(:,2:3)-PT(:,4:5)).^2,2));
        distO = sqrt(sum((O(:,2:3)-O(:,4:5)).^2,2));
        plot(1:N,distP,'-','Color',[1,0,0],'DisplayName','distance Particle');
        plot(1:N,distPT,'-','Color',[0,0,0],'DisplayName','distance True');
        plot(1:N,distO,'-','Color',[0,0,1],'DisplayName','distance observed');
        yl = ylim();
        yl(1)=0;
        stairs((1:N),.98*yl(2)*P(:,1)+0.01,'m-','DisplayName','Predicted State');
        stairs(1:N,.95*yl(2)*PT(:,1)+0.01,'m--','DisplayName','True State');
%         ylim(yl);
        legend('location','best');
        title('Distance/State');

        subplot(2,2,4);
        cla();
        hold('on');
        llh_all = smc.computeLLH_debug(O',P');
        llhT_all = smc.computeLLH_debug(O',PT');
        plot(1:N,llh_all(1,:),'k-','DisplayName','llh z_t');
        plot(1:N,llhT_all(1,:),'k--','DisplayName','llh z_t (true)');
        plot(1:N,sum(llh_all(2:5,:)),'-','Color',[1,0,0],'DisplayName','llh y_t');
        plot(1:N,sum(llhT_all(2:5,:)),'--','Color',[1,0,0],'DisplayName','llh y_t (true)');
        plot(1:N,sum(llh_all(6:9,:)),'-','Color',[0,1,0],'DisplayName','llh x_t');
        plot(1:N,sum(llhT_all(6:9,:)),'--','Color',[0,1,0],'DisplayName','llh x_t (true)');
        legend('location','best');
        title('LLH');
        fprintf('Reported LLH: %.6g \\ %.6g TrueParticle: %.6g\n',PsLLH,sum(llh_all(:)),sum(llhT_all(:)));
      
%         plot(1:N,llh_all(2,:),'-','Color',[1,0,0],'DisplayName','llh ax_t');
%         plot(1:N,llhT_all(2,:),'--','Color',[1,0,0],'DisplayName','llh ax_t (true)');
%         plot(1:N,llh_all(3,:),'-','Color',[.5,0,0],'DisplayName','llh ay_t');
%         plot(1:N,llhT_all(3,:),'--','Color',[.5,0,0],'DisplayName','llh ay_t (true)');
%         plot(1:N,llh_all(4,:),'-','Color',[0,0,1],'DisplayName','llh ax_t');
%         plot(1:N,llhT_all(4,:),'--','Color',[0,0,1],'DisplayName','llh ax_t (true)');
%         plot(1:N,llh_all(5,:),'-','Color',[0,0,.5],'DisplayName','llh ay_t');
%         plot(1:N,llhT_all(5,:),'--','Color',[0,0,.5],'DisplayName','llh ay_t (true)');

    end            
end        

