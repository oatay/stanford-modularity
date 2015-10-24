clear all; close all;
w1 = 1; w2 = 1;
a = 1; b = 1;
k1 = 0.1; k2 = 0.1;
h = 2;
g = 1; k3 = 0.01; w3 = 1;


logarray = zeros(102,1);
logarray(1) = 1e-10;
for i = 2:102
    logarray(i) = logarray(i-1)*1.4;
end

ystar1 = logarray;

xstar1 = w1./(k1^h+ystar1.^h);

xstar2 = logarray;
ystar2 = w2./(k2^h+xstar2.^h);

figure(1)

loglog(xstar1,ystar1,xstar2,ystar2,'LineWidth',2)
hold on
plot(1,1,'o')
plot(1e2,1e-4,'x')
plot(1e-4,1e2,'x')

figure(2)

xdom = linspace(1e-5,1e3,10);
ydom = linspace(1e-5,1e3,10);

[X,Y] = meshgrid(xdom,ydom); % generate mesh of domain

U = w1./(k1*k1+Y.*Y)-X; % dx/dt
V = w2./(k2*k2+X.*X)-Y;   % dy/dt

s_init = [50 20];
[t1,s]= ode45(@(t,y) XYstate(t,y,a,b,k1,k2,w1,w2,h),[0,50],s_init);
plot(t1,s(:,1),t1,s(:,2),'LineWidth',2)

figure(3)
hold on
xlim([-5 110])
ylim([-5 110])
plot(xstar1,ystar1,xstar2,ystar2,'LineWidth',2)


plot(1e2,1e-4,'x')
plot(1e-4,1e2,'x')
u1 = gradient(s(:,1),t1);
v1 = gradient(s(:,2),t1);
quiver(s(:,1),s(:,2),u1,v1,'r')

hold off
% hold on
% quiver(X,Y,U,V)
% 
% hold off

figure(4)
%p_init = [1 0.1 0.8];
p_init = [50 20 10];
%options = odeset('NonNegative', 1:3);
[t,p] = ode45(@(t,y) XYZoscillator(t,y,a,b,g,k1,k2,k3,w1,w2,w3,h),[0,50],p_init);
plot(t,p(:,1),t,p(:,2),t,p(:,3),'LineWidth',2)

figure(5)
x1=1e2; x2 = 1e-4;
semilogy(1:100,x1.*ones(100,1),1:100,x2.*ones(100,1),'LineWidth',2);
ylim([1e-6,1e4])


figure(6)
xdom = linspace(1e-5,1e2,10);
ydom = linspace(1e-5,1e2,10);
zdom = linspace(1e-5,40,10);

[X2,Y2,Z2] = meshgrid(xdom,ydom,zdom); % generate mesh of domain

U2 = -a.*(X2-w2./(k2*k2+Y2.^h));
V2 = -b.*(Y2-(w1./(k1*k1+X2.^h)).*(w3./(k3*k3+Z2.^h)));
T2 = -g.*(Z2-w1./(k1*k1+X2.^h));

quiver3(X2,Y2,Z2,U2,V2,T2,1000)

figure(7)
U3 = gradient(p(:,1),t);
V3 = gradient(p(:,2),t);
Z3 = gradient(p(:,3),t);
scale = 1;
quiver3(p(:,1),p(:,2),p(:,3),U3,V3,Z3,scale)

quiver(p(:,1),p(:,2),U3,V3,scale)