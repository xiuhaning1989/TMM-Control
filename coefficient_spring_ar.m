function u =coefficient_spring_ar(x1,x2,y1,y2,theta,alpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 global psi_bb; 
global a_r; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bb=theta+psi_bb;


ar_c=a_r*cos(t_bb);     ar_s=a_r*sin(t_bb);


%%

u=(sqrt((-x1+x2+ar_c)^2+(y1-y2+ar_s)^2)-a_r)/...
    sqrt((-x1+x2+ar_c)^2+(y1-y2+ar_s)^2);

