%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=Gauss(A,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metodo de Gauss
% Se eligio triangular superior
% para llamar la funcion, guardarla como: Gauss
% para llamarla: Gauss(A,b), Debe definirse A y b.
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
b=transpose(b);
%%%%%%%%%%%%%%%%%%%
n = size(A,1);
%%%%%%%%%%%%%%%%%%%
for k=1:n
    for j=k+1:n
    lambda=A(j,k)/A(k,k);
   A(j,k:n) = A(j,k:n) - lambda*A(k,k);
    b(j)=b(j)-lambda*b(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = n:-1:1
    b(k) = (b(k) - A(k,k+1:n)*b(k+1:n))/A(k,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solucion
x=b
end