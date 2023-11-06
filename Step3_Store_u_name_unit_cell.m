clear all
clc
load Coor_hexagon.mat


%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
% n=5;m=5;

%%
%Store names of variables of entire unit cells: x and y replaced by 1 and 2,respectively, in the vectors
%u_ijx1, u_ijy1, u_ijx2, u_ijy2, u_ijx3, u_ijy3

U_entire_name=[];U_bottom_name=[];U_left_name=[];U_top_name=[];U_right_name=[];

%Store the first row, i=1, only point C:  u_1jx3, u_1jy3
i=1;
for j=1:m
     u=zeros(2,1);
     result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));   u(1) = str2num(result);
     result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));   u(2) = str2num(result);
     if j<10
         u=u/1e3;
     elseif j<100
         u=u/1e4;
     else
         u=u/1e5;
     end
        U_entire_name=[U_entire_name;u];
        U_bottom_name=[U_bottom_name;u]; % store variables on the bottom boundary
end

%Store 2 to n+1 row, point A, B, C: u_ijx1, u_ijy1, u_ijx2, u_ijy2, u_ijx3, u_ijy3,
for i=2:n+1
    for j=1:m
        u=zeros(6,1);
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(1));   u(1) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(1));   u(2) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(2));   u(3) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(2));   u(4) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));   u(5) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));   u(6) = str2num(result);
        if j<10
            u=u/1e3;
        elseif j<100
            u=u/1e4;
        else
            u=u/1e5;
        end
        U_entire_name=[U_entire_name;u];
        if j==1
           U_left_name=[U_left_name;u(3:4)]; 
        end
        if i==n+1
           U_top_name=[U_top_name;u(5:6)];  
        end
    end
end

%Store the last column, j=m+1, from 2 to n+1 row, only point B:  u_1jx2, u_1jy2
j=m+1;
for i=2:n+1
     u=zeros(2,1);
     result = strcat(num2str(i),num2str(j),num2str(1),num2str(2));   u(1) = str2num(result);
     result = strcat(num2str(i),num2str(j),num2str(2),num2str(2));   u(2) = str2num(result);
     if j<10
         u=u/1e3;
     elseif j<100
         u=u/1e4;
     else
         u=u/1e5;
     end
        U_entire_name=[U_entire_name;u];
        U_right_name=[U_right_name;u];   % store variables on the right boundary
end

length(U_entire_name)
6*n*m+2*n+2*m

save('Variables_names_m_by_n_lattice','U_entire_name','U_bottom_name',...
'U_left_name','U_top_name','U_right_name') 
