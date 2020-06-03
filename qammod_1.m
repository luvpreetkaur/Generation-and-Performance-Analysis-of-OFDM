function y = qammod_1(x, m)

k = sqrt(m);
r = 2*(0:k-1) - k + 1;
[xi, yi] = meshgrid(r);
c = xi + 1j*flipud(yi);
y = c(x+1);


