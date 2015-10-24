%clear all; close all; clc;
%THE PARAMETERS USED IN THE PAPER
t0 = 120; tf  = 500; % time-window
aF =3; % given alpha-factor
aF2 = aF; % alpha at the second part
n = 8;
hold all
IC = zeros(8,1);
A = 200; % max. ste12 concentration
C = 1; % alpha factor half-maximum for fus3pp
G = 200; % max. fus3pp concentration
IC(1) = G / (C + exp(-log10(aF))); % fus3pp
IC(2) = 60; % initial total far1
T = 10; % determines half-maximum for transcription
IC(3) = A*(0.35 + 0.65*(1-exp(-aF/T))); % ste12 
IC(4) = 0; % activated far1
IC(5) = 0; % activated cln1/2-clb5/6, >400 after G1/S
IC(6) = 0; % inactivated cln1/2-clb5/6
IC(8) = 30; % initial cell size
IC(7) = IC(8)/2; % Cln3 / 1 nM = ~25 molecules/cell

% parameters
k = zeros(1,17);
k(1) = .01; % Ste12 produces k1*25 Far1 molecules/min
k(2) = .01; % Fus3 activation of Far1 - fast  *******
k(3) = .1; % active Far1 dephosphorylation
k(4) = .05; % active Cln1/2 segregation of Far1-act % strong ****
k(5) = .01; % separation from Cln1/2-Far1 complex
k(6) = 50; % Cln1/2 self-activation production rate % fast
k(7) = 1; % Cln3 activation of Cln1/2 production rate % strong since small Cln3
k(8) = 0.01; % Far1-tot degradation
k(9) = 0.1; % size increase rate
k(10) = .01; % far1act from inactive complex (cln12 degraded)
k(11) = .1; % cln1/2act from inactive complex (far1 degraded)
k(12) = .1; % Cln3 degradation rate

k(13) = 2; % nuclear to cytoplasmic far1 ratio
k(14) = 5; % Cln1/2-Clb5/6 half-maximal activation conc. ******
k(15) = 10; % Cln3 half-maximal activation conc.
k(16) = 100; % half maximal size
k(17) = 10; % the max Cln3 conc. 


[t, output] = ode45(@(t,y) far1clb5mod(t,y,k,n), [t0,tf], IC);

% decrease af 
t0= 0;
IC = output(size(output,1),:);
IC(1) = G / (C + exp(-log10(aF2))); % fus3pp
IC(3) = A*(0.35 + 0.65*(1-exp(-aF2/T))); % ste12 
[t2, output2] = ode45(@(t,y) far1clb5mod(t,y,k,n), [t0, tf], IC);
output = [output; output2];
t2 = t(end) + t2;
t = [t; t2];

figure(1)
subplot(2,4,1)
plot(t,output(:,4))
title('active far1')
subplot(2,4,2)
plot(t,output(:,5))
title('active cln1/2')
subplot(2,4,3)
plot(t,output(:,6))
title('inactive cln12far1')
subplot(2,4,4)
plot(t,output(:,7))
title('cln3 vs. time')
subplot(2,4,5)
plot(t,output(:,8))
title('size')
subplot(2,4,6)
plot(output(:,4),output(:,5))
title('cln12 vs far1')
subplot(2,4,7)
plot(t,output(:,2))
title('total far1')
subplot(2,4,8)
plot(t,output(:,2)-output(:,4))
title('nonphosphorylated far1')

figure(2)
plot(t,output(:,4),'LineWidth',2)
title('Effect of [aF] on time to cell cycle re-entry', 'FontSize',12)
xlabel('time (min)')
ylabel('active Far1 concentration (nM)')
xlim([0 1000])
ylim([0 200])
hold all
tentry = t(find(output(:,5) > 250, 1));
plot([160 tentry],[max(output(:,4)) max(output(:,4))])
plot([tentry tentry],[0 max(output(:,4))])
plot([160 160], [0 max(output(:,4))])

