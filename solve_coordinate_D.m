function [U]=solve_coordinate_D(alpha,theta,gamma)

U=zeros(6,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; global b_b; global c_b;
global psi_ab; global psi_bb; global psi_cb;
global a_r; global b_r; global c_r;
global psi_ar; global psi_br; global psi_cr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U(1)=c_b-c_r*cos(gamma+psi_ab+psi_br)+a_b*cos(gamma-alpha+psi_ab+psi_br)-...
    b_r*cos(theta+psi_bb+psi_cr)+b_b*cos(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)-...
    a_r*cos(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);

U(2)=-c_r*sin(gamma+psi_ab+psi_br)+a_b*sin(gamma-alpha+psi_ab+psi_br)+...
    b_r*sin(theta+psi_bb+psi_cr)-b_b*sin(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)-...
    a_r*sin(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);

U(3)=b_r*sin(theta+psi_bb+psi_cr)-b_b*sin(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)-...
    a_r*sin(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);

U(4)=c_r*sin(gamma+psi_ab+psi_br)-a_b*sin(gamma-alpha+psi_ab+psi_br)+...
    a_r*sin(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);

U(5)=b_r*cos(theta+psi_bb+psi_cr)-b_b*cos(theta+psi_bb+psi_cr+alpha+psi_cb+psi_ar)+...
    a_r*cos(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);

U(6)=-c_r*cos(gamma+psi_ab+psi_br)+a_b*cos(gamma-alpha+psi_ab+psi_br)-...
    a_r*cos(gamma-theta-alpha-psi_bb-psi_cr-psi_cb-psi_ar);
