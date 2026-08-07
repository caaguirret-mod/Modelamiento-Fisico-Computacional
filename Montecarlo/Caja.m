%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p]=Caja(N,Nl)
%Montecarlo en una caja. N=total, NL=particulas a la izq
%Nr=particulas a la derecha
%grafica la cantidad de particulas para t, diferentes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%N=1000; %total
%Nl=900; %inicial izq
Nr=N-Nl; %inicial der 
tMax=N.^2;
x=rand(1,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:N

  p(t)=Nl./N;
   x=rand(1);
if p<x
    Nl = Nl - 1; %Disminuye una
else
    Nr = Nr + 1; % Aumenta una
end 

end 
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1);
plot(p,'b--','LineWidth',2); grid 


end 