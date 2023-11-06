clear all
clc
load Homogeneous_lattice_angles.mat
load Variables_names_m_by_n_lattice.mat
load Coor_hexagon.mat
% load Variables_names_unit_cell.mat

%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
% n and m are from Coor_hexagon.mat

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

global l_s; l_s=(a_b+b_r)/2*0.95;

global k_bond; global k_spring; 
global eta_bond; global eta_spring; 
k_bond=1e4; k_spring=100;
eta_bond=100; eta_spring=20;

global m_b; global m_r;
m_b=0.1;m_r=0.1;

%%

%Choose the angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344)
% i_alpha is from Coor_hexagon.mat
alpha=Alpha(i_alpha);
gamma=Gamma(i_alpha);
theta=Theta(i_alpha);

theta_matrix=theta*ones(n+1,m);
alpha_matrix=alpha*ones(n+1,m);


%%
% Define the y-direction of the bottom edge are always zero. Lattice can
% slide at the bottom edge1
result = strcat(num2str(2),num2str(1),num2str(1),num2str(1));   U_fix_name(1,:) = str2num(result);
result = strcat(num2str(2),num2str(1),num2str(2),num2str(1));   U_fix_name(2,:) = str2num(result);
result = strcat(num2str(2),num2str(1),num2str(1),num2str(2));   U_fix_name(3,:) = str2num(result);
result = strcat(num2str(2),num2str(1),num2str(2),num2str(2));   U_fix_name(4,:) = str2num(result);
result = strcat(num2str(2),num2str(2),num2str(1),num2str(2));   U_fix_name(5,:) = str2num(result);
result = strcat(num2str(2),num2str(2),num2str(2),num2str(2));   U_fix_name(6,:) = str2num(result);
if j<10
    U_fix_name=U_fix_name/1e3;
elseif j<100
    U_fix_name=U_fix_name/1e4;
else
    U_fix_name=U_fix_name/1e5;
end

U_fix_name=[U_bottom_name(1:end);U_fix_name];  

% Define external force at the designated point of the lattice
f_U_name=[];
for j=1:1
result = strcat(num2str(n+1),num2str(j),num2str(2),num2str(3));   f_U_name1 = str2num(result);
if j<10
    f_U_name1=f_U_name1/1e3;
elseif j<100
    f_U_name1=f_U_name1/1e4;
else
    f_U_name1=f_U_name1/1e5;
end
f_U_name=[f_U_name;f_U_name1];
end
F_external=[0].*ones(length(f_U_name),1);


delta_t=1e-4; T=30;
Time=0:delta_t:T;
omega=2*pi;
Force_ext=F_external.*exp(-1*(Time-1).^2).*cos(omega*Time);

u=zeros(2*(6*n*m+2*n+2*m),1);

for i = 1:length(Time)-1
    
    k1 = equation_motion_unit_cell_dissipation(Time(i),u(:,i),n,m,alpha_matrix,theta_matrix,...
        U_entire_name,U_fix_name,Force_ext(:,i),f_U_name);
    k2 = equation_motion_unit_cell_dissipation(Time(i)+0.5*delta_t,u(:,i)+0.5*k1*delta_t,n,m,alpha_matrix,theta_matrix,...
        U_entire_name,U_fix_name,Force_ext(:,i),f_U_name);    
    k3 = equation_motion_unit_cell_dissipation(Time(i)+0.5*delta_t,u(:,i)+0.5*k2*delta_t,n,m,alpha_matrix,theta_matrix,...
        U_entire_name,U_fix_name,Force_ext(:,i),f_U_name); 
    k4 = equation_motion_unit_cell_dissipation(Time(i)+delta_t,u(:,i)+k3*delta_t,n,m,alpha_matrix,theta_matrix,...
        U_entire_name,U_fix_name,Force_ext(:,i),f_U_name); 
    u(:,i+1) = u(:,i)+((k1+2*k2+2*k3+k4)/6)*delta_t;  
end


save('Dynamic_2_1_k2=100_eta2=20_f=0.mat')






















