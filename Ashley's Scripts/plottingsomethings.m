% -- to make heat map of first cargo appearance time
% load Ashley File
load('E:\Zapalog Project\GluA1release_REs\060118\merges\slip1\1_merge_AshleyFile.mat');
h = aa.plot_DHFR_minFrame;
% to change the colorbar
h.Children(5).Label.String = 'hI GUYS';
h.Children(6).XLabel.String = 'some test info';
h.Children(6).XLabel.Visible = 'on'
h.Children(6).Title.String= 'My Title'
h.Children(6).Title.FontSize = 18;
h.Children(6).Title.Visible= 'on '


%%
ff= figure; % make a blank figure to plot into 
hist3data = distVtimemat(distVtimemat(:,2)>=50,:); % this is the data to be plotted;
% hist3data = distVtimemat;
hist3(hist3data,'CDataMode','auto','Nbins',[max(distVtimemat(:,1)) 30],'FaceColor','interp');
xlabel('Time in Frames','FontSize',14);
ylabel('Distance From Soma (microns)','FontSize',14);
