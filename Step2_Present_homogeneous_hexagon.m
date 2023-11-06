clear all
clc
% load Homogeneous_lattice_angles.mat
load Homogeneous_lattice_angles_bistable_new.mat

%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
n=10;m=10;
Coor_hexagon_x=ones(n-1,m-1,7); 
Coor_hexagon_y=ones(n-1,m-1,7);
% x, y coordinates for each hexagon from A, B, C, D, E, F, to A 
Coor_unit_cell_x=ones(n+1,m,9);
Coor_unit_cell_y=ones(n+1,m,9);
% x, y coordinates for each unit cell from C, A, B, C, A', B' C, P, to Q

%Choose the angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344)
i_alpha=16;
alpha=Alpha(i_alpha);
gamma=Gamma(i_alpha);
theta=Theta(i_alpha);

% alpha=Alphanon(i_alpha);
% gamma=Gammanon(i_alpha);
% theta=Thetanon(i_alpha);

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



%%
%Define transformation vectors of coordinates

%From i,j to i+1,j+1
shift_diagonal=[c_b-c_r*cos(gamma+psi_ab+psi_br)+a_b*cos(gamma-alpha+psi_ab+psi_br),...
    -c_r*sin(gamma+psi_ab+psi_br)+a_b*sin(gamma-alpha+psi_ab+psi_br)];

%From i,1, to i+1,1
kappa_vertical=3*pi-alpha-theta-psi_ar-psi_cr-psi_bb;
shift_vertical=[b_r*cos(theta+psi_bb+psi_cr)+a_b*cos(kappa_vertical),...
    -b_r*sin(theta+psi_bb+psi_cr)+a_b*sin(kappa_vertical)];
% phi=gamma+theta+psi_ab+psi_bb-2*pi; %phi=0 for homogeneous lattice

%From 1,j to 1,j+1
kappa_horizontal=psi_ab+gamma-pi;
shift_horizontal=[c_b+a_r*cos(kappa_horizontal),a_r*sin(kappa_horizontal)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Represent an unit cell in a local coordinate system
A = [0,0];
B = [-a_r*cos(psi_ab+gamma-pi),-a_r*sin(psi_ab+gamma-pi)];  
C = [b_r*cos(theta+psi_bb+psi_cr), -b_r*sin(theta+psi_bb+psi_cr)]; 
Aa = [b_r*cos(theta+psi_bb+psi_cr)+a_b*cos(3*pi-alpha-theta-psi_ar-psi_cr-psi_bb),...
    -b_r*sin(theta+psi_bb+psi_cr)+a_b*sin(3*pi-alpha-theta-psi_ar-psi_cr-psi_bb)];
Bb = [b_r*cos(theta+psi_bb+psi_cr)-b_b*cos(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar),...
    -b_r*sin(theta+psi_bb+psi_cr)+b_b*sin(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)];
P= [b_r/2*cos(theta+psi_bb+psi_cr), -b_r/2*sin(theta+psi_bb+psi_cr)];
Q= [b_r*cos(theta+psi_bb+psi_cr)-a_b/2*cos(alpha+theta+psi_ar+psi_cr+psi_bb),...
    -b_r*sin(theta+psi_bb+psi_cr)+a_b/2*sin(alpha+theta+psi_ar+psi_cr+psi_bb)];

Coor_base_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]-shift_vertical;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce n (row) by m (column) unit cells

for j=1:m % generate the 1st row of lattice (bottom edge)
    Coor_unit_cell_x(1,j,:)=Coor_base_unit_cell(:,1)+(j-1)*shift_horizontal(1);  
    Coor_unit_cell_y(1,j,:)=Coor_base_unit_cell(:,2)+(j-1)*shift_horizontal(2);
end

for i=1:n+1 %generate the 1st column of the lattice (left edge)
    Coor_unit_cell_x(i,1,:)=Coor_base_unit_cell(:,1)+(i-1)*shift_vertical(1);  
    Coor_unit_cell_y(i,1,:)=Coor_base_unit_cell(:,2)+(i-1)*shift_vertical(2);
end

for j=1:m %generate unit cells diagonally from the bottom edge
    if m<n+1
        for k=0:m-j
            Coor_unit_cell_x(1+k,j+k,:)=Coor_unit_cell_x(1,j,:)+k*shift_diagonal(1);
            Coor_unit_cell_y(1+k,j+k,:)=Coor_unit_cell_y(1,j,:)+k*shift_diagonal(2); 
        end
    else
        if j<=max(n+1,m)-min(n+1,m)
            for k=0:min(n+1,m)-1
                Coor_unit_cell_x(1+k,j+k,:)=Coor_unit_cell_x(1,j,:)+k*shift_diagonal(1);
                Coor_unit_cell_y(1+k,j+k,:)=Coor_unit_cell_y(1,j,:)+k*shift_diagonal(2); 
            end
        else
            for k=0:max(n+1,m)-j
                Coor_unit_cell_x(1+k,j+k,:)=Coor_unit_cell_x(1,j,:)+k*shift_diagonal(1);
                Coor_unit_cell_y(1+k,j+k,:)=Coor_unit_cell_y(1,j,:)+k*shift_diagonal(2); 
            end
        end
    end
end

for i=1:n+1 %generate unit cells diagonally from the left edge
    if m<n+1
        if i<=max(n+1,m)-min(n+1,m)
            for k=0:min(n+1,m)-1
                Coor_unit_cell_x(i+k,1+k,:)=Coor_unit_cell_x(i,1,:)+k*shift_diagonal(1);
                Coor_unit_cell_y(i+k,1+k,:)=Coor_unit_cell_y(i,1,:)+k*shift_diagonal(2);
            end
        else
            for k=0:max(n+1,m)-i
                Coor_unit_cell_x(i+k,1+k,:)=Coor_unit_cell_x(i,1,:)+k*shift_diagonal(1);
                Coor_unit_cell_y(i+k,1+k,:)=Coor_unit_cell_y(i,1,:)+k*shift_diagonal(2);
            end
        end
    else
        for k=0:n+1-i
            Coor_unit_cell_x(i+k,1+k,:)=Coor_unit_cell_x(i,1,:)+k*shift_diagonal(1);
            Coor_unit_cell_y(i+k,1+k,:)=Coor_unit_cell_y(i,1,:)+k*shift_diagonal(2);
        end
    end
end

%Plot all unit cells to form a n by m lattice
angle=-atan2(shift_horizontal(2),shift_horizontal(1));
rotation_kappa=[cos(angle),-sin(angle);
            sin(angle),cos(-angle)];
        
figure;
for i=1:n+1
    for j=1:m
        for k=1:9
            XY_unit(1,k)=Coor_unit_cell_x(i,j,k);XY_unit(2,k)=Coor_unit_cell_y(i,j,k);
        end
        XY_unit=rotation_kappa*XY_unit;
        if i==1
            hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'b')
            hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'b-','linewidth',0.8);
        elseif i==n+1           
            hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'r')
            hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'r-','linewidth',0.8);
        else
            hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'r')
            hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'r-','linewidth',0.8);
            hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'b')
            hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'b-','linewidth',0.8);
            hold on;plot(XY_unit(1,8:9),XY_unit(2,8:9),'g-','linewidth',1.5)
            axis equal
        end
    end
end
axis equal

% figure;
% for i=1:n+1
%     for j=1:m
%         for k=1:9
%             XY_unit(1,k)=Coor_unit_cell_x(i,j,k);XY_unit(2,k)=Coor_unit_cell_y(i,j,k);
%         end
%         XY_unit=rotation180*XY_unit;
%         if i==1
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
%             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',0.1);
%         elseif i==n+1           
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
%             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',0.1);
%         else
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
%             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',0.1);
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
%             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',0.1);
% %             hold on;plot(XY_unit(1,8:9),XY_unit(2,8:9),'g-','linewidth',2)
%         end
%     end
% end
% axis equal

title(['Homogeneous Lattice of ' num2str(n) ' rows and '...
    num2str(m) ' columns with angle \alpha = ' num2str(alpha)])

save('Coor_hexagon.mat','Coor_unit_cell_x','Coor_unit_cell_y','alpha',...
    'gamma','theta','n','m','i_alpha','rotation_kappa') 

%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Represent nodes of based hexagon in a local coordinate system
% A = [0,0];
% B = [c_b,0];  
% C = [c_b-c_r*cos(gamma+psi_ab+psi_br),-c_r*sin(gamma+psi_ab+psi_br)];
% D=[c_b-c_r*cos(gamma+psi_ab+psi_br)+a_b*cos(gamma-alpha+psi_ab+psi_br),...
%     -c_r*sin(gamma+psi_ab+psi_br)+a_b*sin(gamma-alpha+psi_ab+psi_br)];
% E = [b_r*cos(theta+psi_bb+psi_cr)-b_b*cos(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar),...
%     -b_r*sin(theta+psi_bb+psi_cr)+b_b*sin(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)];
% F = [b_r*cos(theta+psi_bb+psi_cr), -b_r*sin(theta+psi_bb+psi_cr)]; 
% 
% Coor_base_hexagon=[A;B;C;D;E;F;A];  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Reproduce n (row) by m (column) hexagons
% 
% for j=1:m-1 % generate the 1st row of lattice (bottom edge)
%     Coor_hexagon_x(1,j,:)=Coor_base_hexagon(:,1)+(j-1)*shift_horizontal(1);  
%     Coor_hexagon_y(1,j,:)=Coor_base_hexagon(:,2)+(j-1)*shift_horizontal(2);
% end
% 
% for i=1:n-1   %generate the 1st column of the lattice (left edge)
%     Coor_hexagon_x(i,1,:)=Coor_base_hexagon(:,1)+(i-1)*shift_vertical(1);  
%     Coor_hexagon_y(i,1,:)=Coor_base_hexagon(:,2)+(i-1)*shift_vertical(2);
% end
% 
% for j=1:m-1  %generate hexagons diagonally from the bottom edge
%     if m-1<n-1
%         for k=0:m-1-j
%             Coor_hexagon_x(1+k,j+k,:)=Coor_hexagon_x(1,j,:)+k*shift_diagonal(1);  
%             Coor_hexagon_y(1+k,j+k,:)=Coor_hexagon_y(1,j,:)+k*shift_diagonal(2);  
%         end
%     else
%         if j<=max(n-1,m-1)-min(n-1,m-1)
%             for k=0:min(n-1,m-1)-1
%                 Coor_hexagon_x(1+k,j+k,:)=Coor_hexagon_x(1,j,:)+k*shift_diagonal(1);  
%                 Coor_hexagon_y(1+k,j+k,:)=Coor_hexagon_y(1,j,:)+k*shift_diagonal(2);  
%             end
%         else
%             for k=0:max(n-1,m-1)-j
%                 Coor_hexagon_x(1+k,j+k,:)=Coor_hexagon_x(1,j,:)+k*shift_diagonal(1);  
%                 Coor_hexagon_y(1+k,j+k,:)=Coor_hexagon_y(1,j,:)+k*shift_diagonal(2);  
%             end
%         end        
%     end
% end
% 
% for i=1:n-1  %generate hexagons diagonally from the left edge
%     if m-1<n-1
%         if i<=max(n-1,m-1)-min(n-1,m-1)
%             for k=0:min(n-1,m-1)-1
%                 Coor_hexagon_x(i+k,1+k,:)=Coor_hexagon_x(i,1,:)+k*shift_diagonal(1);  
%                 Coor_hexagon_y(i+k,1+k,:)=Coor_hexagon_y(i,1,:)+k*shift_diagonal(2); 
%             end
%         else
%             for k=0:max(n-1,m-1)-i
%                 Coor_hexagon_x(i+k,1+k,:)=Coor_hexagon_x(i,1,:)+k*shift_diagonal(1);  
%                 Coor_hexagon_y(i+k,1+k,:)=Coor_hexagon_y(i,1,:)+k*shift_diagonal(2); 
%             end
%         end
%     else
%         for k=0:n-1-i
%             Coor_hexagon_x(i+k,1+k,:)=Coor_hexagon_x(i,1,:)+k*shift_diagonal(1);  
%             Coor_hexagon_y(i+k,1+k,:)=Coor_hexagon_y(i,1,:)+k*shift_diagonal(2); 
%         end
%     end
% end
% 
% %Plot all unit cells to form a n by m lattice
% figure(2);
% for i=1:n-1
%     for j=1:m-1
%         for k=1:7
%             X_hex(k)=Coor_hexagon_x(i,j,k);Y_hex(k)=Coor_hexagon_y(i,j,k);
%         end
%         hold on; plot(X_hex,Y_hex,'k-o','linewidth',2,'MarkerSize',3)
%     end
% end
% title(['Homogeneous Lattice of ' num2str(n) ' rows and '...
%     num2str(m) ' columns with angle \alpha = ' num2str(alpha)])
% 
% 
% %%
% %Plot all unit cells to form a n by m lattice
% figure(3);
% for i=1:n+1
%     for j=1:m
%         for k=1:9
%             X_unit(k)=Coor_unit_cell_x(i,j,k);Y_unit(k)=Coor_unit_cell_y(i,j,k);
%         end
%         if i==1
%             hold on;fill(X_unit(4:7),Y_unit(4:7),'b')
%             hold on;plot(X_unit(4:7),Y_unit(4:7),'b-','linewidth',0.2);
%         elseif i==n+1           
%             hold on;fill(X_unit(1:4),Y_unit(1:4),'r')
%             hold on;plot(X_unit(1:4),Y_unit(1:4),'r-','linewidth',0.2);
%         else
%             hold on;fill(X_unit(1:4),Y_unit(1:4),'r')
%             hold on;plot(X_unit(1:4),Y_unit(1:4),'r-','linewidth',0.2);
%             hold on;fill(X_unit(4:7),Y_unit(4:7),'b')
%             hold on;plot(X_unit(4:7),Y_unit(4:7),'b-','linewidth',0.2);
%             hold on;plot(X_unit(8:9),Y_unit(8:9),'k-','linewidth',2)
%         end
%     end
% end
% 
% for i=1:n-1
%     for j=1:m-1
%         for k=1:7
%             X_hex(k)=Coor_hexagon_x(i,j,k);Y_hex(k)=Coor_hexagon_y(i,j,k);
%         end
%         hold on; plot(X_hex,Y_hex,'g-o','linewidth',2,'MarkerSize',3)
%     end
% end
% title(['Homogeneous Lattice of ' num2str(n) ' rows and '...
%     num2str(m) ' columns with angle \alpha = ' num2str(alpha)])
