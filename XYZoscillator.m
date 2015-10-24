function dpdt = XYZoscillator(t,p,a,b,g,k1,k2,k3,w1,w2,w3,h)
    dpdt = zeros(3,1);
    dpdt(1) = -a*(p(1)-w2/(k2*k2+p(2)^h));
    dpdt(2) = -b*(p(2)-(w1/(k1*k1+p(1)^h))*(w3/(k3*k3+p(3)^h)));
    dpdt(3) = -g*(p(3)-w1/(k1*k1+p(1)^h));
    
