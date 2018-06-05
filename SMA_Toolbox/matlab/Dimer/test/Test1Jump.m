
close all

%% Setup some parameters

R_real=.05;  %(microns)Seperation of dimers 
S_R=1e-9;     %Sigma discribing actual spread of seperation
dt=1e-3;
S_D=2*dt;      %Sigma for diffusion
S_T=1*dt;      %Signa for rotational diffusion

%Starting condions:

A0=[0 0];
B0=[R_real 0];


%% Move with bound model:

N=100;
clear A1 B1 LL_F LL_BCT
for nn=1:N
    C0=(A0+B0)/2;
    C1=C0+S_D*randn(1,2);
    X=B0-A0;
    T0=atan2(X(2),X(1));
    T1=T0+S_T*randn;
    A1(nn,:)=C1-R_real/2*[cos(T1) sin(T1)];
    B1(nn,:)=C1+R_real/2*[cos(T1) sin(T1)];
    R=norm(B1(nn,:)-A1(nn,:));
    
    %likelihood free
    LL_F(nn)=sum(log(normpdf(A0-A1(nn,:),0,S_D)))+sum(log(normpdf(B0-B1(nn,:),0,S_D)));
    
    %likelihood bound (C,theta)
%     LL_BCT(nn)=sum(log(normpdf(C1-C0,0,S_D)))+log(normpdf(T1-T0,0,S_T))+log(normpdf(R,R_real,S_R))+log(R);
    LL_BCT(nn)=sum(log(normpdf(C1-C0,0,S_D)))+log(normpdf(T1-T0,0,S_T));
    %Note: 'R' is the Jacaobian from below.  Not sure this is correct
    %usage. 
    
end
A1_B=A1;
B1_B=B1;

[CF,BF]=hist(LL_F,100);
[CB,BB]=hist(LL_BCT,BF);
figure
bar(BF,CF,'b')
hold on
bar(BB,CB,'r')
legend('Free Model','Bound Model')
xlabel('Log Likelihood')
ylabel('Counts')
title('Bound->Bound')

LLD=(LL_BCT-LL_F);
BPrcnt=sum(LLD>0)/N;

%% move with free model


clear A1 B1 LL_F LL_BCT
for nn=1:N
    
    A1(nn,:)=A0+S_D*randn(1,2);
    B1(nn,:)=B0+S_D*randn(1,2);
    R=norm(B1(nn,:)-A1(nn,:));
    
    C1=(A1(nn,:)+B1(nn,:))/2;
    X=B1-A1;
    T1=atan2(X(2),X(1));
    
    %likelihood free
    LL_F(nn)=sum(log(normpdf(A0-A1(nn,:),0,S_D)))+sum(log(normpdf(B0-B1(nn,:),0,S_D)));
    
    %likelihood bound (C,theta,r)
%     LL_BCT(nn)=sum(log(normpdf(C1-C0,0,S_D)))+log(normpdf(T1-T0,0,S_T))+log(normpdf(R,R_real,S_R))+log(R);
    LL_BCT(nn)=sum(log(normpdf(C1-C0,0,S_D)))+log(normpdf(T1-T0,0,S_T));
    %Note: 'R' is the Jacaobian from below.  Not sure this is correct
    %usage. 
    
end
A1_F=A1;
B1_F=B1;

[CF,BF]=hist(LL_F,100);
[CB,BB]=hist(LL_BCT,BF);
figure
bar(BF,CF,'b')
hold on
bar(BB,CB,'r')
legend('Free Model','Bound Model')
xlabel('Log Likelihood')
ylabel('Counts')
title('Bound->Free')

LLD=(LL_BCT-LL_F);
BPrcnt=sum(LLD>0)/N


%% Calculate the Jacobian for change of variables. 

syms r Cx Cy Ax Ay Bx By theta

Ax=Cx-r/2*cos(theta);
Ay=Cy-r/2*sin(theta);
Bx=Cx+r/2*cos(theta);
By=Cy+r/2*sin(theta);

J=jacobian([Ax,Bx,Ay,By],[r,theta,Cx,Cy]);
simplify(det(J));


figure();
plot(A1_B(:,1), A1_B(:,2), 'ro','MarkerFaceColor',[1,0,0],'DisplayName','A Bound');
hold();
plot(B1_B(:,1), B1_B(:,2), 'bo','MarkerFaceColor',[0,0,1],'DisplayName','B Bound');
plot(A1_F(:,1), A1_F(:,2), 'o','MarkerSize',5,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[.5,0,0],'DisplayName','A Free');
plot(B1_F(:,1), B1_F(:,2), 'o','MarkerSize',5,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,.5],'DisplayName','B Free');
axis('equal');
legend('location','best');

