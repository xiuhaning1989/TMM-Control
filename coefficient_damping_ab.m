function u =coefficient_damping_ab(x1,x2,y1,y2,theta,alpha,configuration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_b=configuration(1,1);
psi_bb=configuration(1,5);
psi_ar=configuration(2,4);psi_cr=configuration(2,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;
ab_c=a_b*cos(at_arbbcr);    ab_s=a_b*sin(at_arbbcr);    

%%

u=(-x1+x2+ab_c)^2+(y1-y2+ab_s)^2;

