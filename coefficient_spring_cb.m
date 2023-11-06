function u =coefficient_spring_cb(x1,x2,y1,y2,theta,alpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; global b_b; global c_b;
 global psi_bb; global psi_cb;
global psi_ar;  global psi_cr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;
at_arbbcrcb=alpha+theta+psi_ar+psi_bb+psi_cr+psi_cb;
% a_ar=alpha+psi_ar;


ab_c=a_b*cos(at_arbbcr);    ab_s=a_b*sin(at_arbbcr);    
bb_c=b_b*cos(at_arbbcrcb);  bb_s=b_b*sin(at_arbbcrcb);  

%%

u=(sqrt((-x1+x2+ab_c-bb_c)^2+(y1-y2+ab_s-bb_s)^2)-c_b)/...
    sqrt((-x1+x2+ab_c-bb_c)^2+(y1-y2+ab_s-bb_s)^2);

