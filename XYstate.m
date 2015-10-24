function dydt = XYstate(t,s,a,b,k1,k2,w1,w2,h)
    dydt = zeros(2,1);
    dydt(1) = -a*(s(1)-(w2/(k2*k2+s(2)^h)));
    dydt(2) = -b*(s(2)-(w1/(k1*k1+s(1)^h)));