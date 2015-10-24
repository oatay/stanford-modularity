 function dy = far1clb5(t,y,k,c)
dy = zeros(8,1);
%modified far1 phosphorylation by fus3

k1 = k(1); % Ste12 produces k1*25 Far1 molecules/min
k2 = k(2); % Fus3 activation of Far1 - fast
k3 = k(3); % active Far1 dephosphorylation
k4 = k(4); % active Cln1/2 segregation of Far1-act % strong
k5 = k(5); % separation from Cln1/2-Far1 complex
k6 = k(6); % Cln1/2 self-activation production rate % relatively fast
k7 = k(7); % Cln3 activation of Cln1/2 production rate % strong since small Cln3
k8 = k(8); % Far1-tot degradation
k9 = k(9); % size increase rate
k10 = k(10); % far1act from inactive complex (cln12 degraded)
k11 = k(11); % cln1/2act from inactive complex (far1 degraded)
k12 = k(12); % Cln3 degradation rate

kn = k(13); % nuclear ratio
E = k(14); % Cln1/2-Clb5/6 half-maximal activation conc.
F = k(15); % Cln3 half-maximal activation conc.
K = k(16); % half maximal size
B = k(17); % the max Cln3 conc. 
n = c; % ultrasensitivity of Cln1/2-Clb5/6 activation of themselves
m = 1; % ultrasensitivity of Cln3 activation of Cln1/2


dy(1) = 0;  % fus3pp
dy(2) = k1*y(3) - k8*y(2); % total far1
dy(3) = 0; % ste12 
dy(4) = kn*(k2*(y(2)-y(4))*y(1)- k3*y(4) - k4*y(5)*y(4)+ k5*y(6) + k10*y(6)); % activated far1
dy(5) = (k6*y(5)^n)/(E^n+y(5)^n) + (k7*y(7)^m)/(F^m+y(7)^m) - k4*y(4)*y(5)+ k5*y(6) + k11*y(6); % activated cln1/2
dy(6) = k4*y(4)*y(5) - k5*y(6) - k10*y(6) - k11*y(6); % inactivated far1-cln1/2
dy(7) = B*y(8)/(K+y(8))- k12*y(7); % cln3
dy(8) = k9; % size