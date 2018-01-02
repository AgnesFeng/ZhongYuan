function r = ray_doppler(fd, dt, N)
%最大多普勒频移、
% clear all;
% close all;
% fd = 25;
% dt = 1/25e6;
% N = 100000;
T = N*dt-dt; 
t = 0:dt:T;

M = 20;%入射波数量。
c = sqrt(2/M); 
w = 2*pi*fd; 
x = 0; y = 0;
for n = 1:M
    alpha = (2*pi*n-pi+(2*pi*rand-pi))/(4*M);
    ph1 = 2*pi*rand - pi;
    ph2 = 2*pi*rand - pi;
    x = randn(1, N); y = randn(1, N);
    x = x + c*cos(w*t*cos(alpha) + ph1);
    y = y + c*cos(w*t*sin(alpha) + ph2);
end
    r = sqrt(x.^2 + y.^2)/sqrt(2); 
    dopDB = 10*log10(r.^2);
end 