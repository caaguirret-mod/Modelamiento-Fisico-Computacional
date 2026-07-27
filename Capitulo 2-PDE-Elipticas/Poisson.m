%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u]=Poisson(J)
%Esquema explicito Ecuacion eliptica
%Ecuacion uxx+uyy=-f(x,y)
%Esquema=u(i,j)=0.25*(dx*dy*feval(f,x(i),y(j))+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
% para llamar la funcion, guardarla como: Poisson
% para llamarla: Poisson(J), J=malla espacial
%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0.0; x2=1.0; x(1)=0.0; 
y1=0.0; y2=1.0; y(1)=0.0;  
%J=20;
u=zeros(J,J);
dx=(x2-x1)./J;
dy=(y2-y1)./J;
%%%%%%%%%%%%%%%%%%%%%%%%%
f=@(x,y) sin(pi*x).*sin(pi*y);
%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:J 
          x(i+1)=i*dx;
          u(i,1)=feval(f,x(i),1);
          u(i,J)=feval(f,x(i),J);
  end 

   for j=1:J   
     y(j+1)=j*dy;
    u(1,j)=feval(f,1,y(j));
    u(J,j)=feval(f,J,y(j)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:J-1
  for i=2:J-1 
   u(i,j)=0.25*(dx*dy*feval(f,x(i),y(j))+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
  end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%
surf(u);
end