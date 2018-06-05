
% Testing the Distribution Functions used in the HMM 




%% Free Pairs:
L=3; %minimum separation for analysis

N=1000; %number of locations

%space them over a larger area than L
x=10*L*rand(N,1);
y=10*L*rand(N,1);

%calculate all possible distances:
X=repmat(x,[1 N]);
Y=repmat(y,[1 N]);
D=triu(sqrt((X-X').^2+(Y-Y').^2));
d=(D(D<L));
d=d(d>0);

Deltad=0.1;
d_bins=(0:Deltad:L)
[a b]=hist(d,d_bins)

%check
a_norm=a/length(d);
sum(a_norm)

figure
bar(d_bins,a_norm)
hold on
plot(d_bins,2/L^2*d_bins*Deltad,'r','linewidth',3)
xlabel('Distance')
ylabel('Counts')
legend('Simulation Results','Probability Distribution')

%% Bound Pairs

s_R=.1; %registration error
s_C=.0; %crystal structure error
s_L=.2; %localization error

N=5000;

%assume the actual position of particle 1 is at (0,0)
%the actual position of particle 2 is
X2=randn(N,1)*s_R+randn(N,1)*s_C;
Y2=randn(N,1)*s_R+randn(N,1)*s_C;

%the observed positions are:
X1=randn(N,1)*s_L;
Y1=randn(N,1)*s_L;
X2=X2+randn(N,1)*s_L;
Y2=Y2+randn(N,1)*s_L;

d=sqrt((X1-X2).^2 +(Y1-Y2).^2);

Deltad=0.05;
d_bins=(0:Deltad:L);
[a b]=hist(d,d_bins);

%check
a_norm=a/length(d);
sum(a_norm)

s_all2=s_C^2+s_R^2+2*s_L^2;

figure
bar(d_bins,a_norm)
hold on
plot(d_bins,1/s_all2*exp(-d_bins.^2/(2*s_all2)).*d_bins*Deltad,'r','linewidth',3)
xlabel('Distance')
ylabel('Counts')
legend('Simulation Results','Probability Distribution')

%% Free pairs given a prior separation

% S^2 = s1^2+s2^2+2*D*deltaT
% P(r|d)=1/(2*pi*S2)*exp(-r^2/(2*S2))*Integral[ exp(d*r/S2*sin(theta)) d_theta]

%the integral of exp(d*r/S2*sin(theta)) d_theta isn't analytic 
% so lets try to approximate exp(A*sin(theta))

I_theta=inline('r/(2*pi*S2)*.01*sum( exp(-r^2/(2*S2)-d^2/(2*S2)+r*d/S2*sin((0:.01:2*pi))))','r','d','S2');

d=500
r=(0:1:10^3);
S2=10^2+10^2+4*500
i_sum=zeros(length(r),1);
for ii=1:length(r)
   i_sum(ii)=I_theta(r(ii),d,S2);
end
figure
plot(r,(i_sum))
   

%%  check for a specific prior observed distance
% rotate so that y=0 in observed distances

N=10000; %number of locations
d_o=4000; %observed separation
S_L=50; %nm localization accuracy
Diff=500; %nm^2/frame

%observed
x1=zeros(N,1);
y1=zeros(N,1);
x2=zeros(N,1)+d_o;
y2=zeros(N,1);

%actual
x1=x1+S_L*randn(N,1);
y1=y1+S_L*randn(N,1);
x2=x2+S_L*randn(N,1);
y2=y2+S_L*randn(N,1);

%diffusion step
x1=x1+sqrt(2*Diff)*randn(N,1);
y1=y1+sqrt(2*Diff)*randn(N,1);
x2=x2+sqrt(2*Diff)*randn(N,1);
y2=y2+sqrt(2*Diff)*randn(N,1);

%observed final
x1=x1+S_L*randn(N,1);
y1=y1+S_L*randn(N,1);
x2=x2+S_L*randn(N,1);
y2=y2+S_L*randn(N,1);

%calculate distances:
d=sqrt((x2-x1).^2+(y2-y1).^2);

Deltad=10;
d_bins=(0:Deltad:max(d))
[a b]=hist(d,d_bins)

%check
a_norm=a/length(d);
sum(a_norm)

% calc P_d

I_theta=inline('r/(2*pi*S2)*.01*sum( exp(-r^2/(2*S2)-d^2/(2*S2)+r*d/S2*sin((0:.01:2*pi))))','r','d','S2');

r=(0:1:max(d));
S2=S_L^2+S_L^2+S_L^2+S_L^2+2*Diff
i_sum=zeros(length(r),1);
for ii=1:length(r)
   i_sum(ii)=I_theta(r(ii),d_o,S2);
end

figure
bar(d_bins,a_norm/Deltad)
hold on
plot(r,i_sum,'r','linewidth',3)
xlabel('Distance')
ylabel('Counts')
legend('Simulation Results','Probability Distribution')





