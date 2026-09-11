
close all; clear all; clc; 
rho=@(x,y) 2*pi^2*cos(pi*x)*cos(pi*y);
x1=0;x2=1;
y1=0;y2=1;
J1=20;J2=J1;
dx=(x2-x1)/J1;dy=(y2-x1)/J2;h=dx.*dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:dx:x2;
y=y1:dy:y2;nr1=(J1-1)^2;  
nr2=(J2-1)^2;
A=zeros(nr1,nr2);
J2=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nr1
A(i,i)=4;
  if i+1 <= nr1 & mod(i,J2-1) ~= 0
    A(i,i+1)=-1; 
  end
if i+J1-1 <= nr2 
    A(i,i+J2-1)=-1;
end
if i-1 >= 1 & mod(i-1,J2-1) ~= 0
A(i,i-1)=-1; 
 end
if i-(J1-1) >= 1 
    A(i,i-(J1-1))=-1;
  end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:J2-1
    y(j) = j*h;
    for i=1:J1-1
    x(i) = i*h;
[fij]=feval(rho,x(i),y(j)); % evaluate rho(xi,yj) 
vecF((J1-1)*(j-1)+i)=h^2*fij; 
    end 
end

U = A\vecF';

for j=1:J2-1
    for i=1:J2-1
U2d(i,j)=U((J2-1)*(j-1)+i); % change into 2-D array end
    end
end

%figure(1);
%[th,r]=meshgrid(0:dth/4:the,r1:dr/4:r2);
%x=r.*cos(th);
%y=r.*sin(th);
%surf(U2d);
