function u =coefficient_spring_cr(x1,x2,y1,y2,theta,alpha,configuration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi_bb=configuration.psi_bb;
a_r=configuration.a_r;
b_r=configuration.b_r;
c_r=configuration.c_r;
psi_cr=configuration.psi_cr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bb=theta+psi_bb;
t_bbcr=theta+psi_bb+psi_cr;


ar_c=a_r*cos(t_bb);     ar_s=a_r*sin(t_bb);
br_c=b_r*cos(t_bbcr);   br_s=b_r*sin(t_bbcr);
 

%%

u=(sqrt((x1-x2+ar_c-br_c)^2+(-y1+y2+ar_s-br_s)^2)-c_r)/...
    sqrt((x1-x2+ar_c-br_c)^2+(-y1+y2+ar_s-br_s)^2);

