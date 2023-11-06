clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blue Triangle
global a_b; global b_b; global c_b;
a_b=0.5;b_b=0.7;c_b=1;

global psi_ab; global psi_bb; global psi_cb;
psi_ab=acos((a_b^2-b_b^2-c_b^2)/(-2*b_b*c_b));
psi_bb=acos((b_b^2-a_b^2-c_b^2)/(-2*a_b*c_b));
psi_cb=acos((c_b^2-b_b^2-a_b^2)/(-2*b_b*a_b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Red Triangle
global a_r; global b_r; global c_r;
a_r=0.4;b_r=0.8;c_r=1;

global psi_ar; global psi_br; global psi_cr;
psi_ar=acos((a_r^2-b_r^2-c_r^2)/(-2*b_r*c_r));
psi_br=acos((b_r^2-a_r^2-c_r^2)/(-2*a_r*c_r));
psi_cr=acos((c_r^2-b_r^2-a_r^2)/(-2*b_r*a_r));

rest_length_srping=0.95:0.002:1;



%%

for i=1:length(rest_length_srping)
    
l_s=(a_b+b_r)/2*rest_length_srping(i);
cos_alpha0=(l_s^2-(a_b/2)^2-(b_r/2)^2)/(-a_b*b_r/2);
alpha0_pola(i)=acos(cos_alpha0)-psi_ar;
% alpha0_pola(i)=2*pi-acos(cos_alpha0)-psi_ar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0_pola(i)=120/180*pi;gamma0_pola(i)=180/180*pi; % initial guess for theta and gamma
% theta must be less than 2pi-psi_cr-psi_bb, and gamma must be below 2pi-psi_ab-psi_br

%First Shot
R0=[theta0_pola(i),gamma0_pola(i)]';   
F_coorD = solve_coordinate_D(alpha0_pola(i),R0(1),R0(2));
f_D=[F_coorD(1),F_coorD(2)]';
J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
R1=R0-J_f_D\f_D;
% i=0;
%Netow's method
while norm(R1-R0)>10^-9  %&&  rcond(J_f)>10^-10
%     i=i+1
    R0=R1;
    F_coorD = solve_coordinate_D(alpha0_pola(i),R0(1),R0(2));
    f_D=[F_coorD(1),F_coorD(2)]';
    J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
    R1=R0-J_f_D\f_D;
end

theta0_pola(i)=R1(1); gamma0_pola(i)=R1(2); %solve for theta and gamma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end

Alpha=alpha0_pola;
Theta=theta0_pola;
Gamma=gamma0_pola;
%%

figure;plot(Alpha,Theta,'k--','linewidth', 2)
hold on; plot(Alpha,Gamma,'r-','linewidth', 2)
ylabel('\theta and \gamma')
xlabel('\alpha')
legend('\theta','\gamma','Location','east')

save('Homogeneous_lattice_angles_bistable_new.mat','Alpha','Theta','Gamma','rest_length_srping') 




    
    