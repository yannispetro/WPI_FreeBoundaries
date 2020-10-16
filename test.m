x = -5:0.01:5;
d0 = 0.2;

f = x./((x.^2 + d0^2).^(1.5));

plot(x,f); hold on

f2 = sign(x)./(x.^2+d0^2);

plot(x,f2)