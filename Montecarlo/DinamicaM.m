function [v1,v2]=DinamicaM(dt)
%Aplicacion de la dinamica molecular
%Sistema de 2 particulas
%Tiene los dos potenciales, LJ y EL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dt = 0.01;          
pasos = 2000;       
m = 1.0;            
epsilon = 1.0;      
sigma = 1.0;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = [1.0, 3.1];    
v1 = [0.1, -0.1];   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2 = [1.2, 2.0];    
v2 = [-0.1, 0.1];   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist_r1 = zeros(pasos, 2);
hist_r2 = zeros(pasos, 2);

    r_vec = r1 - r2;
    r = norm(r_vec);
    inv_r = sigma / r;
    inv_r6 = inv_r^6;
    inv_r12 = inv_r6^2;
    F_escalar = (24 * epsilon / r) * (2 * inv_r12 - inv_r6);
    %F_escalar =  (inv_r);
    
    r_unitario = r_vec / r;
    F1 = F_escalar * r_unitario; 
    F2 = -F1;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:pasos
   
    hist_r1(t, :) = r1;
    hist_r2(t, :) = r2;
    r1 = r1 + v1 * dt + 0.5 * (F1 / m) * dt^2;
    r2 = r2 + v2 * dt + 0.5 * (F2 / m) * dt^2;
    
    F1_ant = F1;
    F2_ant = F2;

    r_vec = r1 - r2;
    r = norm(r_vec);
    inv_r = sigma / r;
    inv_r6 = inv_r^6;
    inv_r12 = inv_r6^2;
    F_escalar = (24 * epsilon / r) * (2 * inv_r12 - inv_r6);
    %F_escalar =  (inv_r);
    r_unitario = r_vec / r;
    F1 = F_escalar * r_unitario; 
    F2 = -F1;   
    
    v1 = v1 + 0.5 * ((F1_ant + F1) / m) * dt;
    v2 = v2 + 0.5 * ((F2_ant + F2) / m) * dt;
end


figure('Color', 'w');
plot(hist_r1(:,1), hist_r1(:,2), 'b-', 'LineWidth', 1); hold on;
plot(hist_r2(:,1), hist_r2(:,2), 'r-', 'LineWidth', 1);


plot(r1(1), r1(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(r2(1), r2(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);


axis equal; grid on


end 