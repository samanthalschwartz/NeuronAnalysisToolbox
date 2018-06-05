%% Create a 1D diffusive trajectory for a bound or free pair
% S_B defines PDF for distance of particle 2 from 1
% I think this is exact without motion blur effect

D = 0.03; %um^2/s;
dT = 0.001; % s


S_D = sqrt(2*D*dT);    %Diffusion sigma (micron)
S_B = 0.010;        %Distance sigma (micron)
S_E = 0.030;   %Observation error (micron)

Duration = 2; %(seconds)
T = ceil(Duration/dT);     %Number of time frames
cp= round(T/2);     %Location of change point from B->F
L=1;        %Box size (I'm not enforcing)

%Create the trajectories
% cps = unique(round(T*rand(2,1)))';
cps = T/2;
[A,B] = simulateDimerInteractions(T, cps, 1);

% A_True=cumsum(S_D*randn(T,1));
% A = A_True+S_E*randn(T,1);
% B = [A_True(1:cp) + S_B*randn(cp,1)+S_E*randn(cp,1);...  %Bound model
%      A_True(cp)+cumsum(S_D*randn(T-cp,1))+S_E*randn(T-cp,1)];%Free model
%      



%% Probability of Dimer

%1: Free
%2: Bound
%O: Observed trajectory

%P(M2|O) = P(O|M2)/(P(O|M2)+P(O|M1)) = 1/(1+P(O|M1)/P(O|M2))
S_DE = sqrt(S_D^2+S_E^2);
S_BE = sqrt(S_B^2+2*S_E^2);

%Change these variables to work on only a subsequence of the observation
seq_beg=1;
seq_end=T;

%Formulas for P(O|M1) and P(O|M2) and P(O|M1)/P(O|M2) - not directly used
%Pratio = 1/L*prod(normpdf(B(2:end),B(1:end-1),S_DE))/prod(normpdf(A,B,S_BE));
% Pratio = 1/L*exp(sum(log(normpdf(B(seq_beg+1:seq_end),B(seq_beg:seq_end-1),S_DE)))-sum(log(normpdf(A(seq_beg:seq_end),B(seq_beg:seq_end),S_BE))));


% Look at P(M2|O) vs length of interaction
N = seq_end-seq_beg;
P1 = zeros(N,1);
P2 = zeros(N,1);
Pratio = zeros(N,1);
PM2 = zeros(N,1);
for nn=seq_beg+1:seq_end
    P1(nn-seq_beg) = -2*L + sum(log(normpdf(A(seq_beg+1:nn),A(seq_beg:nn-1),S_DE))) + sum( log( normpdf(B(seq_beg+1:nn),B(seq_beg:nn-1),S_DE) ) );
    P2(nn-seq_beg) = -L + sum(log(normpdf(A(seq_beg+1:nn),A(seq_beg:nn-1),S_DE))) + sum(log(normpdf(A(seq_beg:nn),B(seq_beg:nn),S_BE)));
    Pratio(nn-seq_beg) = 1/L*exp(sum(log(normpdf(B(seq_beg+1:nn),B(seq_beg:nn-1),S_DE)))-sum(log(normpdf(A(seq_beg:nn),B(seq_beg:nn),S_BE))));
    PM2(nn-seq_beg) = 1/(1+Pratio(nn-seq_beg));
end

figure();
%Show trajectory
Tvec = (1:T);
subplot(3,1,1);
plot(Tvec,A,'r')
hold on
plot(Tvec,B,'b')

for cp = cps
    plot([cp cp],ylim(),'k:');
end
xlabel('X Position')
ylabel('Time')


%Show the Probability changing over length
ObsVec=(seq_beg+1:seq_end);
subplot(3,1,2);
plot(ObsVec,PM2,'linewidth',3);
hold on
plot(ObsVec,1-PM2,'linewidth',3);
for cp = cps
    plot([cp cp],ylim(),'k:');
end
% plot(ObsVec, P2-P1);
xlabel('Observations')
ylabel('Probability')
legend('Bound Probablilty','Free Probabililty')
ylim([0 1]);

subplot(3,1,3)
plot(ObsVec, P1);
hold on;
plot(ObsVec, P2);
xlabel('Observation');
ylabel('LLH');

%Degub Check code
% figure
% plot(ObsVec,Pratio)
% xlabel('Observations')
% ylabel('P1/P2')





