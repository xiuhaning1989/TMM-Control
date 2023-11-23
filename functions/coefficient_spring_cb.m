function u = coefficient_spring_cb(x1, x2, y1, y2, theta, alpha, configuration)

    a_b = configuration.a_b;
    b_b = configuration.b_b;
    c_b = configuration.c_b;
    psi_bb = configuration.psi_bb;
    psi_cb = configuration.psi_cb;
    psi_ar = configuration.psi_ar;
    psi_cr = configuration.psi_cr;

    at_arbbcr = alpha + theta + psi_ar + psi_bb + psi_cr;
    at_arbbcrcb = alpha + theta + psi_ar + psi_bb + psi_cr + psi_cb;

    ab_c = a_b * cos(at_arbbcr);
    ab_s = a_b * sin(at_arbbcr);
    bb_c = b_b * cos(at_arbbcrcb);
    bb_s = b_b * sin(at_arbbcrcb);

    u = (sqrt((-x1 + x2 + ab_c - bb_c) ^ 2 + (y1 - y2 + ab_s - bb_s) ^ 2) - c_b) / ...
        sqrt((-x1 + x2 + ab_c - bb_c) ^ 2 + (y1 - y2 + ab_s - bb_s) ^ 2);

end
