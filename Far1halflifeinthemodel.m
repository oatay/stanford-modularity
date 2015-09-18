% Far1 half-life modular vs. non-modular

k = 0.01;
tau = 1/k;
halflife = log(2)*tau;
yvals = [halflife halflife 0 0];
xvals = [180 400 401 420];
figure(1)
plot(xvals, yvals)
ylim([0 90])
xlim([150 500])

t = 1:1:500;
k_nonmod = k.*t./50;
tau_nonmod = 1./k_nonmod;
halflife_nonmod = log(2).*tau_nonmod;
figure(2)
plot(t(220:500),halflife_nonmod(21:301))
ylim([0 200])
xlim([180 500])