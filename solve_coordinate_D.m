function [U] = solve_coordinate_D(configuration, alpha, theta, gamma)

    U = zeros(6, 1);

    a_b = configuration.a_b;
    b_b = configuration.b_b;
    c_b = configuration.c_b;
    psi_ab = configuration.psi_ab;
    psi_bb = configuration.psi_bb;
    psi_cb = configuration.psi_cb;
    a_r = configuration.a_r;
    b_r = configuration.b_r;
    c_r = configuration.c_r;
    psi_ar = configuration.psi_ar;
    psi_br = configuration.psi_br;
    psi_cr = configuration.psi_cr;

    U(1) = c_b - c_r * cos(gamma + psi_ab + psi_br) + a_b * cos(gamma - alpha + psi_ab + psi_br) - ...
        b_r * cos(theta + psi_bb + psi_cr) + b_b * cos(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar) - ...
        a_r * cos(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);

    U(2) = -c_r * sin(gamma + psi_ab + psi_br) + a_b * sin(gamma - alpha + psi_ab + psi_br) + ...
        b_r * sin(theta + psi_bb + psi_cr) - b_b * sin(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar) - ...
        a_r * sin(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);

    U(3) = b_r * sin(theta + psi_bb + psi_cr) - b_b * sin(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar) - ...
        a_r * sin(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);

    U(4) = c_r * sin(gamma + psi_ab + psi_br) - a_b * sin(gamma - alpha + psi_ab + psi_br) + ...
        a_r * sin(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);

    U(5) = b_r * cos(theta + psi_bb + psi_cr) - b_b * cos(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar) + ...
        a_r * cos(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);

    U(6) = -c_r * cos(gamma + psi_ab + psi_br) + a_b * cos(gamma - alpha + psi_ab + psi_br) - ...
        a_r * cos(gamma - theta - alpha - psi_bb - psi_cr - psi_cb - psi_ar);
end