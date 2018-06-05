
Dphi = 5e3; %rad^2/s
Dc = 0.250; %um^2/s
dt = 1e-3; %s
rho = 0.050; %um

c_0 = [0,0];
phi_0 = 3*pi/4;

a_0 = c_0 + rho/2*[cos(phi_0), sin(phi_0)];
b_0 = c_0 - rho/2*[cos(phi_0), sin(phi_0)];
ms_0 = 5;
ms_1 = 2;

figure('Position',[0,0,600,800]);
subplot(2,1,1);
hold();

N=1000;
phi_1 = phi_0 + randn(N,1)*sqrt(2*Dphi*dt);
c_1 = repmat(c_0,N,1) + randn(N,2)*sqrt(2*Dc*dt);
a_1 = c_1 + rho/2*[cos(phi_1), sin(phi_1)];
b_1 = c_1 - rho/2*[cos(phi_1), sin(phi_1)];






plot([a_1(:,1)';b_1(:,1)'],[a_1(:,2)';b_1(:,2)'],'-','Color',[0.85,0.85,0.85],'LineWidth',0.5);
plot(c_1(:,1), c_1(:,2),'o','MarkerEdgeColor','none','MarkerFaceColor','g','MarkerSize',ms_1);
plot(a_1(:,1), a_1(:,2),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',ms_1);
plot(b_1(:,1), b_1(:,2),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',ms_1);
axis('equal');
plot([a_0(1),b_0(1)], [a_0(2),b_0(2)], '-', 'LineWidth',2,'Color','k');
plot(a_0(1),a_0(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms_0);
plot(b_0(1),b_0(2),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms_0);
plot(c_0(1),c_0(2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',ms_0);

delta_a = a_1 - repmat(a_0,N,1);
delta_b = b_1 - repmat(b_0,N,1);
delta_c = c_1 - repmat(c_0,N,1);

subplot(2,1,2);
histogram(delta_a(:,1),'DisplayStyle','stairs','EdgeColor',[1,0,0],'Normalization','probability');
histogram(delta_a(:,2),'DisplayStyle','stairs','EdgeColor',[0.5,0,0],'Normalization','probability');
hold();
histogram(delta_b(:,1),'DisplayStyle','stairs','EdgeColor',[0,0,1],'Normalization','probability');
% histogram(delta_b(:,2),'DisplayStyle','stairs','EdgeColor',[0,0,0.5],'Normalization','probability');
histogram(delta_c(:,1),'DisplayStyle','stairs','EdgeColor',[0,1,0],'Normalization','probability');
% histogram(delta_c(:,2),'DisplayStyle','stairs','EdgeColor',[0,0.5,0],'Normalization','probability');


% K=100;
% phi_samp = randn(K,1)*sqrt(2*Dphi*dt);
% a_x = c_0(1) + rho/2 * cos(phi_samp);
% 
% a_xs = linspace(c_1(1,:)-a_0(1), c_1(1,:)+a_0(1),100);
% pa_x = acos(2/rho*(a_xs-c_x))

