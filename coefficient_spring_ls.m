function u =coefficient_spring_ls(x1,x2,y1,y2,theta,alpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; 
global psi_bb; 
 global b_r; 
global psi_ar;  global psi_cr;
 global l_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_bbcr=theta+psi_bb+psi_cr;
at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;


br_c=b_r*cos(t_bbcr);   br_s=b_r*sin(t_bbcr);
ab_c=a_b*cos(at_arbbcr);    ab_s=a_b*sin(at_arbbcr);    

%%

u=(sqrt((-x1+x2+br_c-ab_c)^2+(y1-y2+br_s-ab_s)^2)-2*l_s)/...
    (4*sqrt((-x1+x2+br_c-ab_c)^2+(y1-y2+br_s-ab_s)^2));

