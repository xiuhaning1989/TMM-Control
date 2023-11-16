function du = equation_motion_unit_cell_dissipation(t, u, n, m, alpha, theta, configuration, U_entire_name, U_fix_name, time, f, f_U_name)

    du = zeros(2 * (6 * n * m + 2 * n + 2 * m), 1); % a column vector store all variables

    % lattice configureation
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
    l_s = configuration.l_s;
    k_bond = configuration.k_bond;
    k_spring = configuration.k_spring;
    eta_bond = configuration.eta_bond;
    eta_spring = configuration.eta_spring;
    m_b = configuration.m_b;
    m_r = configuration.m_r;

    k1 = k_bond; k2 = k_spring; eta1 = eta_bond; eta2 = eta_spring;
    alpha = [alpha, alpha(:, m)];
    theta = [theta, theta(:, m)];
    ff = interp1(time, f, t);

    % Find all points that are fixed
    for xx = 1:length(U_fix_name)
        p_0(xx) = find(abs(U_entire_name - U_fix_name(xx)) < 1e-5);
        u(2 * p_0(xx) - 1) = 0; u(2 * p_0(xx)) = 0;
    end

    for i = 1:n + 1

        for j = 1:m + 1
            % Create coefficients in the unit cell of ij
            t_bb = theta(i, j) + psi_bb;
            t_bbcr = theta(i, j) + psi_bb + psi_cr;
            at_arbbcr = alpha(i, j) + theta(i, j) + psi_ar + psi_bb + psi_cr;
            at_arbbcrcb = alpha(i, j) + theta(i, j) + psi_ar + psi_bb + psi_cr + psi_cb;

            ar_c = a_r * cos(t_bb); ar_s = a_r * sin(t_bb);
            br_c = b_r * cos(t_bbcr); br_s = b_r * sin(t_bbcr);
            ab_c = a_b * cos(at_arbbcr); ab_s = a_b * sin(at_arbbcr);
            bb_c = b_b * cos(at_arbbcrcb); bb_s = b_b * sin(at_arbbcrcb);

            % Create all names of variables related to unit cell of i,j
            Uij_entire = zeros(20, 1);
            result = strcat(num2str(i), num2str(j), num2str(1), num2str(1)); Uij_entire(1) = str2double(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(1)); Uij_entire(2) = str2double(result);
            result = strcat(num2str(i), num2str(j), num2str(1), num2str(2)); Uij_entire(3) = str2double(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(2)); Uij_entire(4) = str2double(result);
            result = strcat(num2str(i), num2str(j), num2str(1), num2str(3)); Uij_entire(5) = str2double(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(3)); Uij_entire(6) = str2double(result);
            result = strcat(num2str(i + 1), num2str(j), num2str(1), num2str(1)); Uij_entire(7) = str2double(result);
            result = strcat(num2str(i + 1), num2str(j), num2str(2), num2str(1)); Uij_entire(8) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j), num2str(1), num2str(1)); Uij_entire(9) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j), num2str(2), num2str(1)); Uij_entire(10) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j), num2str(1), num2str(3)); Uij_entire(11) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j), num2str(2), num2str(3)); Uij_entire(12) = str2double(result);
            result = strcat(num2str(i), num2str(j + 1), num2str(1), num2str(2)); Uij_entire(13) = str2double(result);
            result = strcat(num2str(i), num2str(j + 1), num2str(2), num2str(2)); Uij_entire(14) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j - 1), num2str(1), num2str(3)); Uij_entire(15) = str2double(result);
            result = strcat(num2str(i - 1), num2str(j - 1), num2str(2), num2str(3)); Uij_entire(16) = str2double(result);
            result = strcat(num2str(i), num2str(j - 1), num2str(1), num2str(1)); Uij_entire(17) = str2double(result);
            result = strcat(num2str(i), num2str(j - 1), num2str(2), num2str(1)); Uij_entire(18) = str2double(result);
            result = strcat(num2str(i + 1), num2str(j + 1), num2str(1), num2str(2)); Uij_entire(19) = str2double(result);
            result = strcat(num2str(i + 1), num2str(j + 1), num2str(2), num2str(2)); Uij_entire(20) = str2double(result);

            if j < 10
                Uij_entire = Uij_entire / 1e3;

                if j == 9
                    Uij_entire(13) = Uij_entire(13) / 10; Uij_entire(14) = Uij_entire(14) / 10;
                    Uij_entire(19) = Uij_entire(19) / 10; Uij_entire(20) = Uij_entire(20) / 10;
                end

            elseif j < 100
                Uij_entire = Uij_entire / 1e4;

                if j == 10
                    Uij_entire(15) = Uij_entire(15) * 10; Uij_entire(16) = Uij_entire(16) * 10;
                    Uij_entire(17) = Uij_entire(17) * 10; Uij_entire(18) = Uij_entire(18) * 10;
                end

            else
                Uij_entire = Uij_entire / 1e5;
            end

            % Equations of motion of uijx/y3
            if i == 1 && j < m + 1 % Points at bottom
                Uij_bottom = [Uij_entire(5); Uij_entire(6); Uij_entire(7); Uij_entire(8); Uij_entire(19); Uij_entire(20)];

                for xx = 1:length(Uij_bottom)
                    p_x(xx) = find(abs(U_entire_name - Uij_bottom(xx)) < 1e-5);
                end

                uijx3 = u(2 * p_x(1) - 1); uijy3 = u(2 * p_x(2) - 1);
                ui1jx1 = u(2 * p_x(3) - 1); ui1jy1 = u(2 * p_x(4) - 1);
                ui1j1x2 = u(2 * p_x(5) - 1); ui1j1y2 = u(2 * p_x(6) - 1);
                duijx3 = u(2 * p_x(1)); duijy3 = u(2 * p_x(2));
                dui1jx1 = u(2 * p_x(3)); dui1jy1 = u(2 * p_x(4));
                dui1j1x2 = u(2 * p_x(5)); dui1j1y2 = u(2 * p_x(6));

                co_damping_ab = ((uijx3 - ui1jx1 + ab_c) * (duijx3 - dui1jx1) + (uijy3 - ui1jy1 - ab_s) * (duijy3 - dui1jy1)) / ...
                    coefficient_damping_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration);
                co_damping_bb = ((uijx3 - ui1j1x2 + bb_c) * (duijx3 - dui1j1x2) + (uijy3 - ui1j1y2 - bb_s) * (duijy3 - dui1j1y2)) / ...
                    coefficient_damping_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration);
                % dx/dt
                du(2 * p_x(1) - 1) = u(2 * p_x(1));
                % dx2/d2t
                du(2 * p_x(1)) = -3 / (m_b) * (k1 * (uijx3 - ui1jx1 + ab_c) * coefficient_spring_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijx3 - ui1j1x2 + bb_c) * coefficient_spring_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijx3 - ui1jx1 + ab_c) * co_damping_ab + eta1 * (uijx3 - ui1j1x2 + bb_c) * co_damping_bb);
                % dy/dt
                du(2 * p_x(2) - 1) = u(2 * p_x(2));
                % dy2/d2t
                du(2 * p_x(2)) = -3 / (m_b) * (k1 * (uijy3 - ui1jy1 - ab_s) * coefficient_spring_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijy3 - ui1j1y2 - bb_s) * coefficient_spring_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijy3 - ui1jy1 - ab_s) * co_damping_ab + eta1 * (uijy3 - ui1j1y2 - bb_s) * co_damping_bb);                
            elseif i == n + 1 && j < m + 1 % Points at top
                Uij_top = [Uij_entire(5); Uij_entire(6); Uij_entire(1); Uij_entire(2); Uij_entire(3); Uij_entire(4)];

                for xx = 1:length(Uij_top)
                    p_x(xx) = find(abs(U_entire_name - Uij_top(xx)) < 1e-5);
                end

                uijx3 = u(2 * p_x(1) - 1); uijy3 = u(2 * p_x(2) - 1);
                uijx1 = u(2 * p_x(3) - 1); uijy1 = u(2 * p_x(4) - 1);
                uijx2 = u(2 * p_x(5) - 1); uijy2 = u(2 * p_x(6) - 1);
                duijx3 = u(2 * p_x(1)); duijy3 = u(2 * p_x(2));
                duijx1 = u(2 * p_x(3)); duijy1 = u(2 * p_x(4));
                duijx2 = u(2 * p_x(5)); duijy2 = u(2 * p_x(6));

                co_damping_br = ((uijx3 - uijx1 + br_c) * (duijx3 - duijx1) + (uijy3 - uijy1 - br_s) * (duijy3 - duijy1)) / ...
                    coefficient_damping_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration);
                co_damping_cr = ((uijx3 - uijx2 - ar_c + br_c) * (duijx3 - duijx2) + (uijy3 - uijy2 + ar_s - br_s) * (duijy3 - duijy2)) / ...
                    coefficient_damping_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration);
                % dx/dt
                du(2 * p_x(1) - 1) = u(2 * p_x(1));
                % dx2/d2t
                du(2 * p_x(1)) = -3 / (m_r) * (k1 * (uijx3 - uijx1 + br_c) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijx3 - uijx2 - ar_c + br_c) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijx3 - uijx1 + br_c) * co_damping_br + eta1 * (uijx3 - uijx2 - ar_c + br_c) * co_damping_cr);
                % dy/dt
                du(2 * p_x(2) - 1) = u(2 * p_x(2));
                % dy2/d2t
                du(2 * p_x(2)) = -3 / (m_r) * (k1 * (uijy3 - uijy1 - br_s) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijy3 - uijy2 + ar_s - br_s) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijy3 - uijy1 - br_s) * co_damping_br + eta1 * (uijy3 - uijy2 + ar_s - br_s) * co_damping_cr);

            elseif j < m + 1 % All other points of uijx/y3 inside
                Uij_inside = [Uij_entire(5); Uij_entire(6); Uij_entire(1); Uij_entire(2);
                            Uij_entire(3); Uij_entire(4); Uij_entire(7); Uij_entire(8); Uij_entire(19); Uij_entire(20)];

                for xx = 1:length(Uij_inside)
                    p_x(xx) = find(abs(U_entire_name - Uij_inside(xx)) < 1e-5);
                end

                uijx3 = u(2 * p_x(1) - 1); uijy3 = u(2 * p_x(2) - 1);
                uijx1 = u(2 * p_x(3) - 1); uijy1 = u(2 * p_x(4) - 1);
                uijx2 = u(2 * p_x(5) - 1); uijy2 = u(2 * p_x(6) - 1);
                ui1jx1 = u(2 * p_x(7) - 1); ui1jy1 = u(2 * p_x(8) - 1);
                ui1j1x2 = u(2 * p_x(9) - 1); ui1j1y2 = u(2 * p_x(10) - 1);
                duijx3 = u(2 * p_x(1)); duijy3 = u(2 * p_x(2));
                duijx1 = u(2 * p_x(3)); duijy1 = u(2 * p_x(4));
                duijx2 = u(2 * p_x(5)); duijy2 = u(2 * p_x(6));
                dui1jx1 = u(2 * p_x(7)); dui1jy1 = u(2 * p_x(8));
                dui1j1x2 = u(2 * p_x(9)); dui1j1y2 = u(2 * p_x(10));

                co_damping_ab = ((uijx3 - ui1jx1 + ab_c) * (duijx3 - dui1jx1) + (uijy3 - ui1jy1 - ab_s) * (duijy3 - dui1jy1)) / ...
                    coefficient_damping_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration);
                co_damping_bb = ((uijx3 - ui1j1x2 + bb_c) * (duijx3 - dui1j1x2) + (uijy3 - ui1j1y2 - bb_s) * (duijy3 - dui1j1y2)) / ...
                    coefficient_damping_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration);
                co_damping_br = ((uijx3 - uijx1 + br_c) * (duijx3 - duijx1) + (uijy3 - uijy1 - br_s) * (duijy3 - duijy1)) / ...
                    coefficient_damping_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration);
                co_damping_cr = ((uijx3 - uijx2 - ar_c + br_c) * (duijx3 - duijx2) + (uijy3 - uijy2 + ar_s - br_s) * (duijy3 - duijy2)) / ...
                    coefficient_damping_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration);
                % dx/dt
                du(2 * p_x(1) - 1) = u(2 * p_x(1));
                % dx2/d2t
                du(2 * p_x(1)) = -3 / (m_r + m_b) * (k1 * (uijx3 - uijx1 + br_c) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijx3 - uijx2 - ar_c + br_c) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijx3 - ui1jx1 + ab_c) * coefficient_spring_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijx3 - ui1j1x2 + bb_c) * coefficient_spring_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijx3 - uijx1 + br_c) * co_damping_br + eta1 * (uijx3 - uijx2 - ar_c + br_c) * co_damping_cr + ...
                    eta1 * (uijx3 - ui1jx1 + ab_c) * co_damping_ab + eta1 * (uijx3 - ui1j1x2 + bb_c) * co_damping_bb);
                % dy/dt
                du(2 * p_x(2) - 1) = u(2 * p_x(2));
                % dy2/d2t
                du(2 * p_x(2)) = -3 / (m_r + m_b) * (k1 * (uijy3 - uijy1 - br_s) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijy3 - uijy2 + ar_s - br_s) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijy3 - ui1jy1 - ab_s) * coefficient_spring_ab(ui1jx1, uijx3, ui1jy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    k1 * (uijy3 - ui1j1y2 - bb_s) * coefficient_spring_bb(ui1j1x2, uijx3, ui1j1y2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                    eta1 * (uijy3 - uijy1 - br_s) * co_damping_br + eta1 * (uijy3 - uijy2 + ar_s - br_s) * co_damping_cr + ...
                    eta1 * (uijy3 - ui1jy1 - ab_s) * co_damping_ab + eta1 * (uijy3 - ui1j1y2 - bb_s) * co_damping_bb);
            end

            %
            if i > 1 % Only bottom points from u1j, and all other points are from u2j to u(n+1)j

                % Equations of motion of uijx/y2
                if j == 1 % Points at left
                    Uij_left = [Uij_entire(3); Uij_entire(4); Uij_entire(1); Uij_entire(2); Uij_entire(5); Uij_entire(6)];

                    for xx = 1:length(Uij_left)
                        p_x(xx) = find(abs(U_entire_name - Uij_left(xx)) < 1e-5);
                    end

                    uijx2 = u(2 * p_x(1) - 1); uijy2 = u(2 * p_x(2) - 1);
                    uijx1 = u(2 * p_x(3) - 1); uijy1 = u(2 * p_x(4) - 1);
                    uijx3 = u(2 * p_x(5) - 1); uijy3 = u(2 * p_x(6) - 1);
                    duijx2 = u(2 * p_x(1)); duijy2 = u(2 * p_x(2));
                    duijx1 = u(2 * p_x(3)); duijy1 = u(2 * p_x(4));
                    duijx3 = u(2 * p_x(5)); duijy3 = u(2 * p_x(6));

                    co_damping_ar = ((uijx2 - uijx1 + ar_c) * (duijx2 - duijx1) + (uijy2 - uijy1 - ar_s) * (duijy2 - duijy1)) / ...
                        coefficient_damping_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration);
                    co_damping_cr = ((uijx2 - uijx3 + ar_c - br_c) * (duijx2 - duijx3) + (uijy2 - uijy3 - ar_s + br_s) * (duijy2 - duijy3)) / ...
                        coefficient_damping_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration);
                    % dx/dt
                    du(2 * p_x(1) - 1) = u(2 * p_x(1));
                    % dx2/d2t
                    du(2 * p_x(1)) = -3 / (m_r) * (k1 * (uijx2 - uijx1 + ar_c) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijx2 - uijx3 + ar_c - br_c) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijx2 - uijx1 + ar_c) * co_damping_ar + eta1 * (uijx2 - uijx3 + ar_c - br_c) * co_damping_cr);
                    % dy/dt
                    du(2 * p_x(2) - 1) = u(2 * p_x(2));
                    % dy2/d2t
                    du(2 * p_x(2)) = -3 / (m_r) * (k1 * (uijy2 - uijy1 - ar_s) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijy2 - uijy3 - ar_s + br_s) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijy2 - uijy1 - ar_s) * co_damping_ar + eta1 * (uijy2 - uijy3 - ar_s + br_s) * co_damping_cr);

                elseif j == m + 1 % Points at right
                    Uij_right = [Uij_entire(3); Uij_entire(4); Uij_entire(15); Uij_entire(16); Uij_entire(17); Uij_entire(18)];

                    for xx = 1:length(Uij_right)
                        p_x(xx) = find(abs(U_entire_name - Uij_right(xx)) < 1e-5);
                    end

                    uijx2 = u(2 * p_x(1) - 1); uijy2 = u(2 * p_x(2) - 1);
                    ui_1j_1x3 = u(2 * p_x(3) - 1); ui_1j_1y3 = u(2 * p_x(4) - 1);
                    uij_1x1 = u(2 * p_x(5) - 1); uij_1y1 = u(2 * p_x(6) - 1);
                    duijx2 = u(2 * p_x(1)); duijy2 = u(2 * p_x(2));
                    dui_1j_1x3 = u(2 * p_x(3)); dui_1j_1y3 = u(2 * p_x(4));
                    duij_1x1 = u(2 * p_x(5)); duij_1y1 = u(2 * p_x(6));

                    co_damping_bb = ((uijx2 - ui_1j_1x3 - bb_c) * (duijx2 - dui_1j_1x3) + (uijy2 - ui_1j_1y3 + bb_s) * (duijy2 - dui_1j_1y3)) / ...
                        coefficient_damping_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration);
                    co_damping_cb = ((uijx2 - uij_1x1 + ab_c - bb_c) * (duijx2 - duij_1x1) + (uijy2 - uij_1y1 - ab_s + bb_s) * (duijy2 - duij_1y1)) / ...
                        coefficient_damping_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration);
                    % dx/dt
                    du(2 * p_x(1) - 1) = u(2 * p_x(1));
                    % dx2/d2t
                    du(2 * p_x(1)) = -3 / (m_b) * (k1 * (uijx2 - ui_1j_1x3 - bb_c) * coefficient_spring_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijx2 - uij_1x1 + ab_c - bb_c) * coefficient_spring_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijx2 - ui_1j_1x3 - bb_c) * co_damping_bb + eta1 * (uijx2 - uij_1x1 + ab_c - bb_c) * co_damping_cb);
                    % dy/dt
                    du(2 * p_x(2) - 1) = u(2 * p_x(2));
                    % dy2/d2t
                    du(2 * p_x(2)) = -3 / (m_b) * (k1 * (uijy2 - ui_1j_1y3 + bb_s) * coefficient_spring_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijy2 - uij_1y1 - ab_s + bb_s) * coefficient_spring_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijy2 - ui_1j_1y3 + bb_s) * co_damping_bb + eta1 * (uijy2 - uij_1y1 - ab_s + bb_s) * co_damping_cb);
                else
                    Uij_inside = [Uij_entire(3); Uij_entire(4); Uij_entire(1); Uij_entire(2);
                                Uij_entire(5); Uij_entire(6); Uij_entire(15); Uij_entire(16); Uij_entire(17); Uij_entire(18)];

                    for xx = 1:length(Uij_inside)
                        p_x(xx) = find(abs(U_entire_name - Uij_inside(xx)) < 1e-5);
                    end

                    uijx2 = u(2 * p_x(1) - 1); uijy2 = u(2 * p_x(2) - 1);
                    uijx1 = u(2 * p_x(3) - 1); uijy1 = u(2 * p_x(4) - 1);
                    uijx3 = u(2 * p_x(5) - 1); uijy3 = u(2 * p_x(6) - 1);
                    ui_1j_1x3 = u(2 * p_x(7) - 1); ui_1j_1y3 = u(2 * p_x(8) - 1);
                    uij_1x1 = u(2 * p_x(9) - 1); uij_1y1 = u(2 * p_x(10) - 1);
                    duijx2 = u(2 * p_x(1)); duijy2 = u(2 * p_x(2));
                    duijx1 = u(2 * p_x(3)); duijy1 = u(2 * p_x(4));
                    duijx3 = u(2 * p_x(5)); duijy3 = u(2 * p_x(6));
                    dui_1j_1x3 = u(2 * p_x(7)); dui_1j_1y3 = u(2 * p_x(8));
                    duij_1x1 = u(2 * p_x(9)); duij_1y1 = u(2 * p_x(10));

                    co_damping_ar = ((uijx2 - uijx1 + ar_c) * (duijx2 - duijx1) + (uijy2 - uijy1 - ar_s) * (duijy2 - duijy1)) / ...
                        coefficient_damping_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration);
                    co_damping_cr = ((uijx2 - uijx3 + ar_c - br_c) * (duijx2 - duijx3) + (uijy2 - uijy3 - ar_s + br_s) * (duijy2 - duijy3)) / ...
                        coefficient_damping_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration);
                    co_damping_bb = ((uijx2 - ui_1j_1x3 - bb_c) * (duijx2 - dui_1j_1x3) + (uijy2 - ui_1j_1y3 + bb_s) * (duijy2 - dui_1j_1y3)) / ...
                        coefficient_damping_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration);
                    co_damping_cb = ((uijx2 - uij_1x1 + ab_c - bb_c) * (duijx2 - duij_1x1) + (uijy2 - uij_1y1 - ab_s + bb_s) * (duijy2 - duij_1y1)) / ...
                        coefficient_damping_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration);
                    % dx/dt
                    du(2 * p_x(1) - 1) = u(2 * p_x(1));
                    % dx2/d2t
                    du(2 * p_x(1)) = -3 / (m_r + m_b) * (k1 * (uijx2 - uijx1 + ar_c) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijx2 - uijx3 + ar_c - br_c) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijx2 - ui_1j_1x3 - bb_c) * coefficient_spring_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijx2 - uij_1x1 + ab_c - bb_c) * coefficient_spring_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijx2 - uijx1 + ar_c) * co_damping_ar + eta1 * (uijx2 - uijx3 + ar_c - br_c) * co_damping_cr + ...
                        eta1 * (uijx2 - ui_1j_1x3 - bb_c) * co_damping_bb + eta1 * (uijx2 - uij_1x1 + ab_c - bb_c) * co_damping_cb);
                    % dy/dt
                    du(2 * p_x(2) - 1) = u(2 * p_x(2));
                    % dy2/d2t
                    du(2 * p_x(2)) = -3 / (m_r + m_b) * (k1 * (uijy2 - uijy1 - ar_s) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijy2 - uijy3 - ar_s + br_s) * coefficient_spring_cr(uijx2, uijx3, uijy2, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijy2 - ui_1j_1y3 + bb_s) * coefficient_spring_bb(uijx2, ui_1j_1x3, uijy2, ui_1j_1y3, theta(i, j), alpha(i, j), configuration) + ...
                        k1 * (uijy2 - uij_1y1 - ab_s + bb_s) * coefficient_spring_cb(uij_1x1, uijx2, uij_1y1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                        eta1 * (uijy2 - uijy1 - ar_s) * co_damping_ar + eta1 * (uijy2 - uijy3 - ar_s + br_s) * co_damping_cr + ...
                        eta1 * (uijy2 - ui_1j_1y3 + bb_s) * co_damping_bb + eta1 * (uijy2 - uij_1y1 - ab_s + bb_s) * co_damping_cb);
                end

                % All points of uijx/y1 are from u21 to u(n+1)m
                if j < m + 1 
                    %Equations of motion of uijx/y1
                    if i == 2 % Points at bottom
                        Uij_inside = [Uij_entire(1); Uij_entire(2); Uij_entire(3); Uij_entire(4); Uij_entire(5); Uij_entire(6);
                                    Uij_entire(7); Uij_entire(8); Uij_entire(11); Uij_entire(12); Uij_entire(13); Uij_entire(14); ];

                        for xx = 1:length(Uij_inside)
                            p_x(xx) = find(abs(U_entire_name - Uij_inside(xx)) < 1e-5);
                        end

                        uijx1 = u(2 * p_x(1) - 1); uijy1 = u(2 * p_x(2) - 1);
                        uijx2 = u(2 * p_x(3) - 1); uijy2 = u(2 * p_x(4) - 1);
                        uijx3 = u(2 * p_x(5) - 1); uijy3 = u(2 * p_x(6) - 1);
                        ui1jx1 = u(2 * p_x(7) - 1); ui1jy1 = u(2 * p_x(8) - 1);
                        ui_1jx3 = u(2 * p_x(9) - 1); ui_1jy3 = u(2 * p_x(10) - 1);
                        uij1x2 = u(2 * p_x(11) - 1); uij1y2 = u(2 * p_x(12) - 1);
                        duijx1 = u(2 * p_x(1)); duijy1 = u(2 * p_x(2));
                        duijx2 = u(2 * p_x(3)); duijy2 = u(2 * p_x(4));
                        duijx3 = u(2 * p_x(5)); duijy3 = u(2 * p_x(6));
                        dui1jx1 = u(2 * p_x(7)); dui1jy1 = u(2 * p_x(8));
                        dui_1jx3 = u(2 * p_x(9)); dui_1jy3 = u(2 * p_x(10));
                        duij1x2 = u(2 * p_x(11)); duij1y2 = u(2 * p_x(12));

                        co_damping_ls1 = ((uijx1 - ui1jx1 - br_c + ab_c) * (duijx1 - dui1jx1) + (uijy1 - ui1jy1 + br_s - ab_s) * (duijy1 - dui1jy1)) / ...
                            coefficient_damping_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration);
                        co_damping_ar = ((uijx1 - uijx2 - ar_c) * (duijx1 - duijx2) + (uijy1 - uijy2 + ar_s) * (duijy1 - duijy2)) / ...
                            coefficient_damping_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration);
                        co_damping_br = ((uijx1 - uijx3 - br_c) * (duijx1 - duijx3) + (uijy1 - uijy3 + br_s) * (duijy1 - duijy3)) / ...
                            coefficient_damping_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_ab = ((uijx1 - ui_1jx3 - ab_c) * (duijx1 - dui_1jx3) + (uijy1 - ui_1jy3 + ab_s) * (duijy1 - dui_1jy3)) / ...
                            coefficient_damping_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_cb = ((uijx1 - uij1x2 - ab_c + bb_c) * (duijx1 - duij1x2) + (uijy1 - uij1y2 + ab_s - bb_s) * (duijy1 - duij1y2)) / ...
                            coefficient_damping_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration);
                        % dx/dt
                        du(2 * p_x(1) - 1) = u(2 * p_x(1));
                        % dx2/d2t
                        du(2 * p_x(1)) = -3 / (m_r + m_b) * (k2 * (uijx1 - ui1jx1 - br_c + ab_c) * coefficient_spring_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx2 - ar_c) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx3 - br_c) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - ui_1jx3 - ab_c) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uij1x2 - ab_c + bb_c) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijx1 - ui1jx1 - br_c + ab_c) * co_damping_ls1 + eta1 * (uijx1 - uijx2 - ar_c) * co_damping_ar + eta1 * (uijx1 - uijx3 - br_c) * co_damping_br + ...
                            eta1 * (uijx1 - ui_1jx3 - ab_c) * co_damping_ab + eta1 * (uijx1 - uij1x2 - ab_c + bb_c) * co_damping_cb);
                        % dy/dt
                        du(2 * p_x(2) - 1) = u(2 * p_x(2));
                        % dy2/d2t
                        du(2 * p_x(2)) = -3 / (m_r + m_b) * (k2 * (uijy1 - ui1jy1 + br_s - ab_s) * coefficient_spring_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy2 + ar_s) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy3 + br_s) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - ui_1jy3 + ab_s) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uij1y2 + ab_s - bb_s) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijy1 - ui1jy1 + br_s - ab_s) * co_damping_ls1 + eta1 * (uijy1 - uijy2 + ar_s) * co_damping_ar + eta1 * (uijy1 - uijy3 + br_s) * co_damping_br + ...
                            eta1 * (uijy1 - ui_1jy3 + ab_s) * co_damping_ab + eta1 * (uijy1 - uij1y2 + ab_s - bb_s) * co_damping_cb);

                    elseif i == n + 1 % Points at top
                        Uij_inside = [Uij_entire(1); Uij_entire(2); Uij_entire(3); Uij_entire(4); Uij_entire(5); Uij_entire(6);
                                    Uij_entire(9); Uij_entire(10); Uij_entire(11); Uij_entire(12); Uij_entire(13); Uij_entire(14); ];

                        for xx = 1:length(Uij_inside)
                            p_x(xx) = find(abs(U_entire_name - Uij_inside(xx)) < 1e-5);
                        end

                        uijx1 = u(2 * p_x(1) - 1); uijy1 = u(2 * p_x(2) - 1);
                        uijx2 = u(2 * p_x(3) - 1); uijy2 = u(2 * p_x(4) - 1);
                        uijx3 = u(2 * p_x(5) - 1); uijy3 = u(2 * p_x(6) - 1);
                        ui_1jx1 = u(2 * p_x(7) - 1); ui_1jy1 = u(2 * p_x(8) - 1);
                        ui_1jx3 = u(2 * p_x(9) - 1); ui_1jy3 = u(2 * p_x(10) - 1);
                        uij1x2 = u(2 * p_x(11) - 1); uij1y2 = u(2 * p_x(12) - 1);
                        duijx1 = u(2 * p_x(1)); duijy1 = u(2 * p_x(2));
                        duijx2 = u(2 * p_x(3)); duijy2 = u(2 * p_x(4));
                        duijx3 = u(2 * p_x(5)); duijy3 = u(2 * p_x(6));
                        dui_1jx1 = u(2 * p_x(7)); dui_1jy1 = u(2 * p_x(8));
                        dui_1jx3 = u(2 * p_x(9)); dui_1jy3 = u(2 * p_x(10));
                        duij1x2 = u(2 * p_x(11)); duij1y2 = u(2 * p_x(12));

                        co_damping_ls2 = ((uijx1 - ui_1jx1 + br_c - ab_c) * (duijx1 - dui_1jx1) + (uijy1 - ui_1jy1 - br_s + ab_s) * (duijy1 - dui_1jy1)) / ...
                            coefficient_damping_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration);
                        co_damping_ar = ((uijx1 - uijx2 - ar_c) * (duijx1 - duijx2) + (uijy1 - uijy2 + ar_s) * (duijy1 - duijy2)) / ...
                            coefficient_damping_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration);
                        co_damping_br = ((uijx1 - uijx3 - br_c) * (duijx1 - duijx3) + (uijy1 - uijy3 + br_s) * (duijy1 - duijy3)) / ...
                            coefficient_damping_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_ab = ((uijx1 - ui_1jx3 - ab_c) * (duijx1 - dui_1jx3) + (uijy1 - ui_1jy3 + ab_s) * (duijy1 - dui_1jy3)) / ...
                            coefficient_damping_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_cb = ((uijx1 - uij1x2 - ab_c + bb_c) * (duijx1 - duij1x2) + (uijy1 - uij1y2 + ab_s - bb_s) * (duijy1 - duij1y2)) / ...
                            coefficient_damping_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration);
                        % dx/dt
                        du(2 * p_x(1) - 1) = u(2 * p_x(1));
                        % dx2/d2t
                        du(2 * p_x(1)) = -3 / (m_r + m_b) * (k2 * (uijx1 - ui_1jx1 + br_c - ab_c) * coefficient_spring_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx2 - ar_c) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx3 - br_c) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - ui_1jx3 - ab_c) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uij1x2 - ab_c + bb_c) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijx1 - ui_1jx1 + br_c - ab_c) * co_damping_ls2 + eta1 * (uijx1 - uijx2 - ar_c) * co_damping_ar + eta1 * (uijx1 - uijx3 - br_c) * co_damping_br + ...
                            eta1 * (uijx1 - ui_1jx3 - ab_c) * co_damping_ab + eta1 * (uijx1 - uij1x2 - ab_c + bb_c) * co_damping_cb);
                        % dy/dt
                        du(2 * p_x(2) - 1) = u(2 * p_x(2));
                        % dy2/d2t
                        du(2 * p_x(2)) = -3 / (m_r + m_b) * (k2 * (uijy1 - ui_1jy1 - br_s + ab_s) * coefficient_spring_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy2 + ar_s) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy3 + br_s) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - ui_1jy3 + ab_s) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uij1y2 + ab_s - bb_s) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijy1 - ui_1jy1 - br_s + ab_s) * co_damping_ls2 + eta1 * (uijy1 - uijy2 + ar_s) * co_damping_ar + eta1 * (uijy1 - uijy3 + br_s) * co_damping_br + ...
                            eta1 * (uijy1 - ui_1jy3 + ab_s) * co_damping_ab + eta1 * (uijy1 - uij1y2 + ab_s - bb_s) * co_damping_cb);
                    else
                        Uij_inside = [Uij_entire(1); Uij_entire(2); Uij_entire(3); Uij_entire(4); Uij_entire(5); Uij_entire(6);
                                    Uij_entire(7); Uij_entire(8); Uij_entire(9); Uij_entire(10); Uij_entire(11); Uij_entire(12);
                                    Uij_entire(13); Uij_entire(14); ];

                        for xx = 1:length(Uij_inside)
                            p_x(xx) = find(abs(U_entire_name - Uij_inside(xx)) < 1e-5);
                        end

                        uijx1 = u(2 * p_x(1) - 1); uijy1 = u(2 * p_x(2) - 1);
                        uijx2 = u(2 * p_x(3) - 1); uijy2 = u(2 * p_x(4) - 1);
                        uijx3 = u(2 * p_x(5) - 1); uijy3 = u(2 * p_x(6) - 1);
                        ui1jx1 = u(2 * p_x(7) - 1); ui1jy1 = u(2 * p_x(8) - 1);
                        ui_1jx1 = u(2 * p_x(9) - 1); ui_1jy1 = u(2 * p_x(10) - 1);
                        ui_1jx3 = u(2 * p_x(11) - 1); ui_1jy3 = u(2 * p_x(12) - 1);
                        uij1x2 = u(2 * p_x(13) - 1); uij1y2 = u(2 * p_x(14) - 1);
                        duijx1 = u(2 * p_x(1)); duijy1 = u(2 * p_x(2));
                        duijx2 = u(2 * p_x(3)); duijy2 = u(2 * p_x(4));
                        duijx3 = u(2 * p_x(5)); duijy3 = u(2 * p_x(6));
                        dui1jx1 = u(2 * p_x(7)); dui1jy1 = u(2 * p_x(8));
                        dui_1jx1 = u(2 * p_x(9)); dui_1jy1 = u(2 * p_x(10));
                        dui_1jx3 = u(2 * p_x(11)); dui_1jy3 = u(2 * p_x(12));
                        duij1x2 = u(2 * p_x(13)); duij1y2 = u(2 * p_x(14));

                        co_damping_ls1 = ((uijx1 - ui1jx1 - br_c + ab_c) * (duijx1 - dui1jx1) + (uijy1 - ui1jy1 + br_s - ab_s) * (duijy1 - dui1jy1)) / ...
                            coefficient_damping_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration);
                        co_damping_ls2 = ((uijx1 - ui_1jx1 + br_c - ab_c) * (duijx1 - dui_1jx1) + (uijy1 - ui_1jy1 - br_s + ab_s) * (duijy1 - dui_1jy1)) / ...
                            coefficient_damping_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration);
                        co_damping_ar = ((uijx1 - uijx2 - ar_c) * (duijx1 - duijx2) + (uijy1 - uijy2 + ar_s) * (duijy1 - duijy2)) / ...
                            coefficient_damping_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration);
                        co_damping_br = ((uijx1 - uijx3 - br_c) * (duijx1 - duijx3) + (uijy1 - uijy3 + br_s) * (duijy1 - duijy3)) / ...
                            coefficient_damping_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_ab = ((uijx1 - ui_1jx3 - ab_c) * (duijx1 - dui_1jx3) + (uijy1 - ui_1jy3 + ab_s) * (duijy1 - dui_1jy3)) / ...
                            coefficient_damping_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration);
                        co_damping_cb = ((uijx1 - uij1x2 - ab_c + bb_c) * (duijx1 - duij1x2) + (uijy1 - uij1y2 + ab_s - bb_s) * (duijy1 - duij1y2)) / ...
                            coefficient_damping_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration);
                        % dx/dt
                        du(2 * p_x(1) - 1) = u(2 * p_x(1));
                        % dx2/d2t
                        du(2 * p_x(1)) = -3 / (m_r + m_b) * (k2 * (uijx1 - ui1jx1 - br_c + ab_c) * coefficient_spring_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration) + ...
                            k2 * (uijx1 - ui_1jx1 + br_c - ab_c) * coefficient_spring_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx2 - ar_c) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uijx3 - br_c) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - ui_1jx3 - ab_c) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijx1 - uij1x2 - ab_c + bb_c) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijx1 - ui1jx1 - br_c + ab_c) * co_damping_ls1 + eta2 * (uijx1 - ui_1jx1 + br_c - ab_c) * co_damping_ls2 + ...
                            eta1 * (uijx1 - uijx2 - ar_c) * co_damping_ar + eta1 * (uijx1 - uijx3 - br_c) * co_damping_br + ...
                            eta1 * (uijx1 - ui_1jx3 - ab_c) * co_damping_ab + eta1 * (uijx1 - uij1x2 - ab_c + bb_c) * co_damping_cb);
                        % dy/dt
                        du(2 * p_x(2) - 1) = u(2 * p_x(2));
                        % dy2/d2t
                        du(2 * p_x(2)) = -3 / (m_r + m_b) * (k2 * (uijy1 - ui1jy1 + br_s - ab_s) * coefficient_spring_ls(uijx1, ui1jx1, uijy1, ui1jy1, theta(i, j), alpha(i, j), configuration) + ...
                            k2 * (uijy1 - ui_1jy1 - br_s + ab_s) * coefficient_spring_ls(ui_1jx1, uijx1, ui_1jy1, uijy1, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy2 + ar_s) * coefficient_spring_ar(uijx1, uijx2, uijy1, uijy2, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uijy3 + br_s) * coefficient_spring_br(uijx1, uijx3, uijy1, uijy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - ui_1jy3 + ab_s) * coefficient_spring_ab(uijx1, ui_1jx3, uijy1, ui_1jy3, theta(i, j), alpha(i, j), configuration) + ...
                            k1 * (uijy1 - uij1y2 + ab_s - bb_s) * coefficient_spring_cb(uijx1, uij1x2, uijy1, uij1y2, theta(i, j), alpha(i, j), configuration) + ...
                            eta2 * (uijy1 - ui1jy1 + br_s - ab_s) * co_damping_ls1 + eta2 * (uijy1 - ui_1jy1 - br_s + ab_s) * co_damping_ls2 + ...
                            eta1 * (uijy1 - uijy2 + ar_s) * co_damping_ar + eta1 * (uijy1 - uijy3 + br_s) * co_damping_br + ...
                            eta1 * (uijy1 - ui_1jy3 + ab_s) * co_damping_ab + eta1 * (uijy1 - uij1y2 + ab_s - bb_s) * co_damping_cb);
                    end
                end % All points of uijx/y1 are from u21 to u(n+1)m

            end % Only bottom points from u1j, and all other points are from u2j to u(n+1)j

        end

    end

    % Acceleration are set to be zero for those points fixed
    % du(2*p_0(xx)-1)=0;
    du(2 * p_0) = 0;

    % Find all points that are applied external forces
    for xx = 1:length(f_U_name)
        p_f(xx) = find(abs(U_entire_name - f_U_name(xx)) < 1e-5);
        du(2 * p_f(xx)) = du(2 * p_f(xx)) + ff(xx);
    end
end