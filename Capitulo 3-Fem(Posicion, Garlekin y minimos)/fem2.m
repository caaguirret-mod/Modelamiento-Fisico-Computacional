%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h]=fem2(Ne)
%Solucion uxx=f(x) via metodo de Garlekin
%u=ax(1-x),como funcion test
%Guardar como Fem2, Ne=numero de elementos que se desean
%llamar fem2(N); en pantalla grafica y presenta h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1.0;          
%Ne = 6;  
Nnod = Ne + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(1)=0;
for k=1:Ne
x(k+1)=x(k)+(1/Ne); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = x(2) - x(1);  
K = zeros(Nnod, Nnod);
F = zeros(Nnod, 1);
f_val = 1.0; 
%%%%%%%%%%%%%%%%%%%%%%%  
for k = 1:Ne
    n1 = k;
    n2 = k + 1;
    K_local = (1./h).*[1, -1; -1, 1];
    F_local = (f_val * h / 2) * [1; 1];
    K(n1:n2, n1:n2) = K(n1:n2, n1:n2) + K_local;
    F(n1:n2) = F(n1:n2) + F_local;
end 
%%%%%%%%%%%%%%%%%%%%%%%
nl = 2:Nnod-1;
u= K(nl,nl) \ F(nl);
u2 = zeros(Nnod, 1);
u(nl) = u;
%%%%%%%%%%%%%%%%%%%%%%%
u(1)=0;
u(Nnod)=0;
%%%%%%%%%%%%%%%%%%%%%%%
plot(x, u, '-o', 'LineWidth', 2);
grid on;
end 