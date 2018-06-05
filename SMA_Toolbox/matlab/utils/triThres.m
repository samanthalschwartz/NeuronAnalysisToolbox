% Mark J. Olah (mjo@cs.unm.edu)
% 07/31/2014

function [Xthresh, threshVal] = triThres(seq, debugPlot)
    %
    % Given an sorted intensities of things that may be signal or background,
    % find a Xthreshold that divides the signal from the background.  We use the
    % triangle method where we find the x value where the distance from the data
    % curve to the triangle enclosing the shape is maximized.
    %
    % See: G. W. Zack et. al. J Histochem Cytochem (25):741 1977
    %
    % [in] seq - sorted sequence of intensity values
    % [in] debugPlot - logical 1=make some plots to help out
    % [out] Xthresh - index of last 'signal' value, all subsequent values are
    % background
    if nargin==1
        debugPlot=0;
    end
    seq=seq(:);
    N=length(seq);
    xs=(0:N-1)';
    seq=double(seq);
    Nseq=seq*(N/seq(1)); %Normalize the sequence so ymax=xmax
    %From (x,seq(x)), we intesect the triangle at (xstar,ystar)
    xstar=(Nseq(1)-Nseq+xs)/2;
    ystar=Nseq(1)-xstar;
    Dstar=sqrt((xs-xstar).^2+(Nseq-ystar).^2); %dist from (x,seq(x)) to (xstar,ystar)
    [~,Xthresh]=max(Dstar); %index of threshold is maximum distance
    threshVal=seq(Xthresh);
    if ~isempty(debugPlot) && (debugPlot>0 || (debugPlot~=0 && ishandle(debugPlot)))
        if ~ishandle(debugPlot)
            figure();
        else
            axes(debugPlot)
        end
        H=Nseq(1)-xs;
        Hp_max=xs-Xthresh+Nseq(Xthresh);
        plot(xs,Nseq,'b-','LineWidth',2.0,'DisplayName','Sorted Data Values');
        hold on
        plot(xs,H,'k-','LineWidth',2.0,'DisplayName','Triangle');
        plot([Xthresh,xstar(Xthresh)],[Hp_max(Xthresh), Hp_max(round(xstar(Xthresh)))],'r-','LineWidth',2.0,...
                'DisplayName','Perp. Line @ Max Dist');
        plot(xs,Dstar,'m-','DisplayName','Dist');
        plot(Xthresh, Dstar(Xthresh),'p','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k','DisplayName','Max Distance');
        plot([Xthresh Xthresh],[0 Nseq(1)],'k:','DisplayName',sprintf('Threshold: %g',threshVal));
        hold off
        legend('Location','East');
        xlabel('Index (n)');
        ylabel('Intensity');
        title('Triangle Method of Background Thresholding');
    end
end

