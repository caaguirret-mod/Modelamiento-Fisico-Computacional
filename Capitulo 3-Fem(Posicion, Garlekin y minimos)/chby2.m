function [h]=chby2(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solucion Txx=-Q;
%metodo espectral de Chebyshev
%Guardala como chby2, se llama como: chby2(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=-1;
x1=-1/2;
x2=1/2;
x3=1;

%N=20;
x(1)=-1;
for j=1:N
x(j+1)=-cos((j*pi)/N);
end 
x0=-1;
x1=-1/2;
x2=1/2;
x3=1;

L0=(x-x1).*(x-x2).*(x-x3);
L0=L0./((x0-x1).*(x0-x2).*(x0-x3));


L1=(x-x0).*(x-x2).*(x-x3);
L1=L1./((x1-x0).*(x1-x2).*(x1-x3));

L2=(x-x0).*(x-x1).*(x-x3);
L2=L2./((x2-x0).*(x2-x1).*(x2-x3));


L3=(x-x0).*(x-x1).*(x-x2);
L3=L3./((x3-x0).*(x3-x1).*(x3-x2));
subplot (2,2,1);
plot(x,L0,'-o', 'LineWidth', 2); hold on; plot(x,L1,'-o', 'LineWidth', 2); hold on;  plot(x,L3,'-o', 'LineWidth', 2);grid
tl=1;
Qx=-1;
Qx2=-2;
tr=1;
dL0dx =@(x)  ((x-x2).*(x-x3) ...
       + (x-x1).*(x-x3) ...
       + (x-x1).*(x-x2)) ...
       ./ ((x0-x1).*(x0-x2).*(x0-x3));

dL1dx =@(x) ((x-x2).*(x-x3) ...
       + (x-x0).*(x-x3) ...
       + (x-x0).*(x-x2)) ...
       ./ ((x1-x0).*(x1-x2).*(x1-x3));

dL2dx =@(x)((x-x1).*(x-x3) ...
       + (x-x0).*(x-x3) ...
       + (x-x0).*(x-x1)) ...
       ./ ((x2-x0).*(x2-x1).*(x2-x3));

dL3dx =@(x) ((x-x1).*(x-x2) ...
       + (x-x0).*(x-x2) ...
       + (x-x0).*(x-x1)) ...
       ./ ((x3-x0).*(x3-x1).*(x3-x2));


D=zeros(4,4);
D=[ feval(dL0dx,x0), feval(dL1dx,x0), feval(dL2dx,x0), feval(dL3dx,x0);
    feval(dL0dx,x1), feval(dL1dx,x1), feval(dL2dx,x1), feval(dL3dx,x1);
    feval(dL0dx,x2), feval(dL1dx,x2), feval(dL2dx,x2), feval(dL3dx,x2);
    feval(dL0dx,x3), feval(dL1dx,x3), feval(dL2dx,x3), feval(dL3dx,x3)]; %Matriz derivada 

D2=D*D; %Matriz segunda derivada
D2(1,1)=1;D2(1,2)=0; D2(1,3)=0; D2(1,4)=0;
D(4,1)=0;D(4,2)=0;D(4,3)=0; D(4,4)=1;
D2(4,:)=0; D2(4,4)=1;
b=[tl;-Qx;-Qx2;tr];
x=D2\b;
subplot(2,2,2);
plot(x,'-o', 'LineWidth', 2); grid
end 