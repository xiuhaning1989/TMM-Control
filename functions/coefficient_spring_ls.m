function u = coefficient_spring_ls(x1, x2, y1, y2, theta, alpha, configuration)

    l_s = configuration.l_s;
    a_b = configuration.a_b;
    b_r = configuration.b_r;
    psi_bb = configuration.psi_bb;
    psi_ar = configuration.psi_ar;
    psi_cr = configuration.psi_cr;

    t_bbcr = theta + psi_bb + psi_cr;
    at_arbbcr = alpha + theta + psi_ar + psi_bb + psi_cr;

    br_c = b_r * cos(t_bbcr);
    br_s = b_r * sin(t_bbcr);
    ab_c = a_b * cos(at_arbbcr);
    ab_s = a_b * sin(at_arbbcr);

    u = (sqrt((-x1 + x2 + br_c - ab_c) ^ 2 + (y1 - y2 + br_s - ab_s) ^ 2) - 2 * l_s) / ...
        (4 * sqrt((-x1 + x2 + br_c - ab_c) ^ 2 + (y1 - y2 + br_s - ab_s) ^ 2));

end
