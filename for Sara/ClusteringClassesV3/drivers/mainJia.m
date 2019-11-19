nSim=25;% how many time to run the simulation
for ii=1:nSim
    [results, results_o, results_c, results_t,results_simulated] = testSLJia();
    N_Collapsed{ii}=results.numobjs_collapsed_pass1;
    N_Col2True{ii}=results_simulated.N_observations;%observations per true localization
end

clear results results_c results_o results_simulated results_t ii
%close all

%over lap the 2 figure
N_Coll = cell2mat(N_Collapsed);
N_C2T  = cell2mat(N_Col2True);
n_Coll = length(N_Coll);
n_C2T  = length(N_C2T);
if n_Coll < n_C2T
   N_Coll = [N_Coll, ones(1, n_C2T - n_Coll)*NaN];
elseif n_Coll > n_C2T
   N_C2T  = [N_C2T,  ones(1, n_Coll - n_C2T)*NaN];
end
N_C = [N_Coll', N_C2T'];
figure()
ah1=axes;
hist(N_C, 40);
%hist(N_Coll, 40);% cyan dots vs. blue dots
hold on
%hist(N_C2T,  40);% cyan dots vs. blue dots
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');%cyan vs blue, red bar
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k');%cyan vs red, blue bar
box off
title({'Simulations'})
ylabel({'Frequency'})
%xlabel({'Collapsed observations per localization'})
legend('localizations','simulated localizations')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(ah1,'LineWidth',3,'FontSize',18)
hold off

%over lap the 2 figure
figure()
ah1=axes;
hist(cell2mat(N_Collapsed), 40);% cyan dots vs. blue dots
hold on
hist(cell2mat(N_Col2True),  40);% cyan dots vs. blue dots
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');%cyan vs blue, red bar
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k');%cyan vs red, blue bar
box off
title({'Simulations'})
ylabel({'Frequency'})
%xlabel({'Collapsed observations per localization'})
legend('localizations','simulated localizations')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(ah1,'LineWidth',3,'FontSize',18)
hold off
