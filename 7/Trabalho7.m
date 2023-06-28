close all
clear all
clc

%% parte I
clear all

dt = 0.01;   %diminuir dt para reduzir aliasing
N = 2^15;   %aumentar N para melhorar resolução (valores + próximos)
t = 0:dt:(N-1)*dt;

domega=2*pi/(N*dt);
omega = -N/2*domega:domega:(N/2-1)*domega;

% Y = exp(-1i*10*t)+exp(1i*20*t);
% Y = sin(t);
Y = fftshift(fft(Y));
y_freq = (dt*abs(Y)).^2;


plot(omega,y_freq,'.-')

%% parte II
clear all

max_z = 4;
max_x = 10;
dz = 0.1;
dx = 0.1;
z = 0:dz:max_z;
x = -max_x:dx:max_x;
N = length(x);

dk = 2*pi/(N*dx);
k = -N/2*dk:dk:(N/2-1)*dk;

Nz = length(z);

q = zeros(Nz, N);


q0 = exp(-x.^2/2);
%q0 = sech(x);


q(1,:) = q0;

q0tr = dk*fftshift(fft(q0));

for i = 2:Nz
    q_tr = q0tr.*exp(-1i/2.*k.^2.*z(i));
    q(i,:) = ifft(ifftshift(q_tr))/dk;
end

mesh(x,z,abs(q).^2)


