function u =coefficient_damping_cb(x1,x2,y1,y2,theta,alpha,configuration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_b=configuration(1,1);b_b=configuration(1,2);
psi_bb=configuration(1,5);psi_cb=configuration(1,6);
psi_ar=configuration(2,4);psi_cr=configuration(2,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;
at_arbbcrcb=alpha+theta+psi_ar+psi_bb+psi_cr+psi_cb;
% a_ar=alpha+psi_ar;


ab_c=a_b*cos(at_arbbcr);    ab_s=a_b*sin(at_arbbcr);    
bb_c=b_b*cos(at_arbbcrcb);  bb_s=b_b*sin(at_arbbcrcb);  

%%

u=(-x1+x2+ab_c-bb_c)^2+(y1-y2+ab_s-bb_s)^2;

