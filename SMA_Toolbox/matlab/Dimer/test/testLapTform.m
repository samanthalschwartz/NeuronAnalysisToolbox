function testLapTform(r,D,t)
    % r - radius
    % D - diffusion constant
    % t - time
    r0s=linspace(0,3*r,100);
    LFS1=@(r0) @(s) (2*r/D)*r0.*besselk(0,r*sqrt(s./D)).*besseli(0,r0.*sqrt(s./D)); %Satomi Eq 12
    LFS2=@(r0) @(s) (2*r/D)*r0.*besselk(0,r0.*sqrt(s./D)).*besseli(0,r*sqrt(s./D)); %Satomi Eq 12 with inverted r,r0
    LFS3=@(r0) @(s) (r/D)*r0.*(besselk(0,r0.*sqrt(s./D)).*besseli(0,r*sqrt(s./D))+besselk(0,r*sqrt(s./D)).*besseli(0,r0.*sqrt(s./D))); %Satomi Eq 13
    LFM=@(r0) @(s) (2*r/D)*r0.*besselk(0,max(r,r0).*sqrt(s./D)).*besseli(0,min(r,r0).*sqrt(s./D)); %Mark Eq. 26&27
    f=@(r0) (r/(D*t)).*r0.*exp(-1/(4*D*t)*(r^2+r0.^2)).*besseli(0,(r/(2*D*t)).*r0);
    fvalsT = arrayfun(f, r0s);
    fvalsS1 = arrayfun(@(r0) laplaceInverse(LFS1(r0),t), r0s);
    fvalsS2 = arrayfun(@(r0) laplaceInverse(LFS2(r0),t), r0s); 
    fvalsS3 = arrayfun(@(r0) laplaceInverse(LFS3(r0),t), r0s); 
    fvalsM = arrayfun(@(r0) laplaceInverse(LFM(r0),t), r0s);
    
    figure();
    plot(r0s,fvalsT,'ko','DisplayName','f(r,r0,t) using direct formula (Satomi 11)');
    hold();
    plot(r0s,fvalsS3,'r-','LineWidth',2,'DisplayName','f(r,r0,t) using (Satomi 13)');
    plot(r0s,fvalsS1,'g-','LineWidth',2,'DisplayName','f(r,r0,t) using (Satomi 12)');
    plot(r0s,fvalsS2,'b-','LineWidth',2,'DisplayName','f(r,r0,t) using (Satomi 12inv)');
    plot(r0s,fvalsM,'m.','MarkerSize',10,'DisplayName','f(r,r0,t) using (Mark Eq.26&27)');
    mv=max([fvalsS1,fvalsS2,fvalsS3,fvalsM]);
    yl=[-mv,mv]*1.1;
    ylim(yl);
    plot([r,r],yl,'k:','DisplayName','r=r0');
    ylabel('f(r,r0,t)');
    xlabel('r0');
    legend('Location','best');
end

