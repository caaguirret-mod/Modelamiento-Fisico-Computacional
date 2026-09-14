close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
k0 = -2;
L = 20;         
Nx = 50;         
x = linspace(-L/2, L/2, Nx)';
dx = x(2) - x(1);
dt=1.0;
Nt=100;
r=dt./(dx^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = -10;         
sigma = 4;
V0 = 10;  
%V=rand(Nx,1);
V = randi(10, Nx, 1);

%for i=1:Nx
%V(i,1)=rand(Nx);
%end
psi = exp(-((x-x0).^2) / (2*sigma^2)) .* exp(-1i * k0 * x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %plot(x,abs(psi))
u=zeros(Nx,Nt);

u(2:Nx-1,1)=psi(2:Nx-1)-1i*r*(psi(3:Nx)+2*psi(2:Nx-1)+psi(1:Nx-2)) +dt.*V(2:Nx-1).*psi(2:Nx-1);

figure('Color', 'w');
for n=2:Nt    
           u(2:Nx-1,n)=u(2:Nx-1,n-1)-1i*r*(u(3:Nx,n-1)+2*u(2:Nx-1,n-1)+u(1:Nx-2,n-1))+...
               +dt.*V(2:Nx-1,1).*u(2:Nx-1,n);
    u(1, :) =0 ; 
    u(end, :) = 0;
 u(1:Nx,n) = u(1:Nx,n)./sqrt(trapz(x, abs(u(1:Nx,n).^2)));

end
% Plot the absolute value of the wave function at the final time step
subplot(2,2,1); 
plot(x, abs(u(:, 30))); grid 

subplot(2,2,2); 
plot(x, abs(u(:, 50))); grid 
subplot(2,2,3); 
plot(x, abs(u(:, 70))); grid 

subplot(2,2,4); 
plot(x, abs(u(:, 100))); grid 
%hPlot = plot(x, abs(u).^2, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
%hold on;
%axis([min(x) max(x) 0 max(abs(psi).^2)*1.2]);
%set(hPlot, 'YData', abs(psi).^2);
%drawnow;