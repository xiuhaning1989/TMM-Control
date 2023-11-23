function u = coefficient_spring_ab(x1, x2, y1, y2, theta, alpha, configuration)
    
    a_b = configuration.a_b;
    psi_bb = configuration.psi_bb;
    psi_ar = configuration.psi_ar;
    psi_cr = configuration.psi_cr;

    at_arbbcr = alpha + theta + psi_ar + psi_bb + psi_cr;
    ab_c = a_b * cos(at_arbbcr); ab_s = a_b * sin(at_arbbcr);

    u = (sqrt((-x1 + x2 + ab_c) ^ 2 + (y1 - y2 + ab_s) ^ 2) - a_b) / ...
        sqrt((-x1 + x2 + ab_c) ^ 2 + (y1 - y2 + ab_s) ^ 2);

end
