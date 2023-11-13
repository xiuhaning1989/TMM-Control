function u =coefficient_spring_ab(x1,x2,y1,y2,theta,alpha,configuration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_b=configuration(1,1);
psi_bb=configuration(1,5);
psi_ar=configuration(2,4);psi_cr=configuration(2,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;

% a_ar=alpha+psi_ar;


ab_c=a_b*cos(at_arbbcr);    ab_s=a_b*sin(at_arbbcr);    

%%

u=(sqrt((-x1+x2+ab_c)^2+(y1-y2+ab_s)^2)-a_b)/...
    sqrt((-x1+x2+ab_c)^2+(y1-y2+ab_s)^2);

