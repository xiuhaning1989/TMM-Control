function u =coefficient_spring_br(x1,x2,y1,y2,theta,alpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 global psi_bb; 
global b_r; 
global psi_cr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_bbcr=theta+psi_bb+psi_cr;




br_c=b_r*cos(t_bbcr);   br_s=b_r*sin(t_bbcr);


%%

u=(sqrt((-x1+x2+br_c)^2+(y1-y2+br_s)^2)-b_r)/...
   sqrt((-x1+x2+br_c)^2+(y1-y2+br_s)^2);

