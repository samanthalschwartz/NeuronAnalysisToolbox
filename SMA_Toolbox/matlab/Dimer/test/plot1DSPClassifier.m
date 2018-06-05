function plot1DSPClassifier(x0, t,  v1, v2)
    % x0 - start position (um)
    % t - time (s)
    % v1 - variance (um^2/s)
    % v2 - variance (um^2/s) [v1 < v2]
    assert(v1<v2);

    eta = @(v1,v2) (v1*v2)/(v2-v1) * log(v2/v1);

    E1 = @(v1,v2) normcdf(-sqrt(eta(v1,v2)/v1),0,1);
    E2 = @(v1,v2) .5 - normcdf(-sqrt(eta(v1,v2)),0,v2);
    E = @(v1,v2) .5 + normcdf(-sqrt(eta(v1,v2)),0,v1) - normcdf(-sqrt(eta(v1,v2)),0,v2);

    figure();
    xs = linspace(x0-5*sqrt(v2), x0+5*sqrt(v2),100);

    ys=normpdf(xs,x0,sqrt(t*v1));
    H=fill([xs,xs(end),xs(1)],[ys,0,0],[1,0,1],'DisplayName','$\mathcal{N}(x_t \mid x_0, tv_1)$');
    H.EdgeColor=[1,0,1];
    H.FaceAlpha=0.4;
    hold();

    ys=normpdf(xs,x0,sqrt(t*v2));
    H=fill([xs,xs(end),xs(1)],[ys,0,0],[0,1,1],'DisplayName','$\mathcal{N}(x_t \mid x_0, tv_2)$');
    H.EdgeColor=[0,1,1];
    H.FaceAlpha=0.4;



    
    xs = linspace(x0-5*sqrt(v2), x0-sqrt(t*eta(v1,v2)),100);
    ys=.5*normpdf(xs,x0,sqrt(t*v1));
    H=fill([xs,xs(end),xs(1)],[ys,0,0],[1,0,0],'DisplayName','Error');
    H.EdgeColor='None';
    H.FaceAlpha=1;

    xs = linspace(x0+sqrt(t*eta(v1,v2)), x0+5*sqrt(v2),100);
    ys=.5*normpdf(xs,x0,sqrt(t*v1));
    H=fill([xs,xs(end),xs(1)],[ys,0,0],[1,0,0],'DisplayName','Error');
    H.EdgeColor='None';
    H.FaceAlpha=1;
    H.Annotation.LegendInformation.IconDisplayStyle='off';

    xs = linspace(x0-sqrt(t*eta(v1,v2)), x0+sqrt(t*eta(v1,v2)),100);
    ys=.5*normpdf(xs,x0,sqrt(t*v2));
    H=fill([xs,xs(end),xs(1)],[ys,0,0],[1,0,0],'DisplayName','Error');
    H.EdgeColor='None';
    H.FaceAlpha=1;
    H.Annotation.LegendInformation.IconDisplayStyle='off';

    yl = ylim();
    plot((x0-sqrt(t*eta(v1,v2))).*[1 1],[0, yl(2)],'k-.','DisplayName','$x_0-\sqrt{\eta t}$');
    plot((x0+sqrt(t*eta(v1,v2))).*[1 1],[0, yl(2)],'k--','DisplayName','$x_0+\sqrt{\eta t}$');

    H=legend('location','best');
    H.Interpreter='latex';
    xlabel('x');
    ylabel('probability');


    figure();
    v2= linspace(1,10,100);
    v2=v2(2:end);
    v10=1;
%     plot(v2,arrayfun(@(v2) eta(1,v2),v2),'DisplayName','$\eta(1,v_2)$');
    hold();
    E1v2=arrayfun(@(v2) E1(v10,v2),v2);
    plot(v2,E1v2,'DisplayName','$E_1(1,v_2)$');
    E2v2=arrayfun(@(v2) E2(v10,v2),v2);
    plot(v2,E2v2,'DisplayName','$E_2(1,v_2)$');
    plot(v2,E1v2+E2v2,'DisplayName','$E(1,v_2)=E_1(1,v_2)+E_2(1,v_2)$');
    Ev2=arrayfun(@(v2) E(v10,v2),v2);
    plot(v2,Ev2,'DisplayName','$E(1,v_2)$');
    xl=xlim();
    yl=ylim();
    ylim([0,yl(2)]);
    plot(xl,[.5, .5],':k','DisplayName','Maximum classification Error = .5')
    H=legend('location','best');
    H.Interpreter='latex';
    xlabel('$v_2$','Interpreter','latex');
end
