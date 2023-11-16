function u = coefficient_damping_bb(x1, x2, y1, y2, theta, alpha, configuration)

    b_b = configuration.b_b;
    psi_bb = configuration.psi_bb;
    psi_cb = configuration.psi_cb;
    psi_ar = configuration.psi_ar;
    psi_cr = configuration.psi_cr;

    at_arbbcrcb = alpha + theta + psi_ar + psi_bb + psi_cr + psi_cb;
    bb_c = b_b * cos(at_arbbcrcb);
    bb_s = b_b * sin(at_arbbcrcb);

    u = (-x1 + x2 + bb_c) ^ 2 + (y1 - y2 + bb_s) ^ 2;

end
