%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=tridag(J)
%Metodo tridiagonal, para solucionar ODE, segundo orden
%P,q,r son funciones de x, x=[0,1];
%frontera y(a)=0.5  y(b)=2.0;
%guardarla como tridag y llamarla tridag(J)
%J maneja la cantidad de pasos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=@(x) 2*pi^2*cos(pi*x);%funcion q
p=@(x) 2*pi^2*sin(pi*x); %funcion p
r=@(x) 2*pi^2*cos(pi*x); %funcion r
x1=0;x2=1;
%J=10;
h=(x2-x1)/J;
h2=h^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:h:x2;
nr1=(J);  
A=zeros(nr1,nr1);
y(1)=-0.5; %valor inicial
y(2)=2.0; %valor inicial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nr1
A(i,i)=2+h2*feval(q,x(i));%diagonales

  if i+1 <= nr1 & mod(i,J) ~= 0
    A(i,i+1)=-1+0.5*h2*feval(p,x(i));%diagonal superior 
  end
if i+J-1 <= nr1 
    A(i,i+J-1)=0;
end
if i-1 >= 1 & mod(i-1,J) ~= 0
A(i,i-1)=-1+0.5*h2*feval(p,x(i+1)); %diagonal inferior 
 end
if i-(J) >= 1 
    A(i,i-(J-1))=0;
  end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(nr1, 1); %Todos menos inicial y final
for i = 2:nr1-1
    b(i) = -h2*feval(r,x(i)); 
end
b(1)=-h2*feval(r,x(1))+(1+0.5*h*feval(p,x(1)))*y(1);
b(nr1)=-h2*feval(r,x(nr1))+(1-0.5*h*feval(p,x(nr1)))*y(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=A\b;
plot(x,'b--','LineWidth',2); grid 
end