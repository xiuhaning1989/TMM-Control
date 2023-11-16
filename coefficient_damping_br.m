function u = coefficient_damping_br(x1, x2, y1, y2, theta, alpha, configuration)

    psi_bb = configuration.psi_bb;
    b_r = configuration.b_r;
    psi_cr = configuration.psi_cr;

    t_bbcr = theta + psi_bb + psi_cr;
    br_c = b_r * cos(t_bbcr);
    br_s = b_r * sin(t_bbcr);

    u = (-x1 + x2 + br_c) ^ 2 + (y1 - y2 + br_s) ^ 2;

end
