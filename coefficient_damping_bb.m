function u =coefficient_damping_bb(x1,x2,y1,y2,theta,alpha,configuration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_b=configuration(1,2);
psi_bb=configuration(1,5);psi_cb=configuration(1,6);
psi_ar=configuration(2,4);psi_cr=configuration(2,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at_arbbcrcb=alpha+theta+psi_ar+psi_bb+psi_cr+psi_cb;
bb_c=b_b*cos(at_arbbcrcb);  bb_s=b_b*sin(at_arbbcrcb);  

%%

u=(-x1+x2+bb_c)^2+(y1-y2+bb_s)^2;

