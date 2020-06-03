function z = qamdemod_1(y, m)
k = sqrt(m);
r = 2*(0:k-1) - k + 1;
[xi, yi] = meshgrid(r);
c = xi + 1j*flipud(yi);
c = c(:);

z = zeros(size(y));

for k = 1:length(y)
    [~, ind] = min(abs(y(k) - c));
    z(k) = ind - 1;
end


