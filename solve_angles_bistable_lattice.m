function [Alpha, Gamma, Theta] = solve_angles_bistable_lattice(configuration)

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

    rest_length_srping = 0.95:0.002:1;

    %%

    for i = 1:length(rest_length_srping)

        l_s = (a_b + b_r) / 2 * rest_length_srping(i);
        cos_alpha0 = (l_s ^ 2 - (a_b / 2) ^ 2 - (b_r / 2) ^ 2) / (-a_b * b_r / 2);
        alpha0_pola(i) = acos(cos_alpha0) - psi_ar;
        % alpha0_pola(i)=2*pi-acos(cos_alpha0)-psi_ar;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta0_pola(i) = 120/180 * pi; gamma0_pola(i) = 180/180 * pi; % initial guess for theta and gamma
        % theta must be less than 2pi-psi_cr-psi_bb, and gamma must be below 2pi-psi_ab-psi_br

        %First Shot
        R0 = [theta0_pola(i), gamma0_pola(i)]';
        F_coorD = solve_coordinate_D(configuration,alpha0_pola(i), R0(1), R0(2));
        f_D = [F_coorD(1), F_coorD(2)]';
        J_f_D = [F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
        R1 = R0 - J_f_D \ f_D;
        % i=0;
        %Netow's method
        while norm(R1 - R0) > 10 ^ -9 %&&  rcond(J_f)>10^-10
            %     i=i+1
            R0 = R1;
            F_coorD = solve_coordinate_D(configuration,alpha0_pola(i), R0(1), R0(2));
            f_D = [F_coorD(1), F_coorD(2)]';
            J_f_D = [F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
            R1 = R0 - J_f_D \ f_D;
        end

        theta0_pola(i) = R1(1); gamma0_pola(i) = R1(2); %solve for theta and gamma
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    Alpha = alpha0_pola;
    Theta = theta0_pola;
    Gamma = gamma0_pola;
    %%

    % figure; plot(Alpha, Theta, 'k--', 'linewidth', 2)
    % hold on; plot(Alpha, Gamma, 'r-', 'linewidth', 2)
    % ylabel('\theta and \gamma')
    % xlabel('\alpha')
    % legend('\theta', '\gamma', 'Location', 'east')

    % save('Homogeneous_lattice_angles_bistable_new.mat', 'Alpha', 'Theta', 'Gamma', 'rest_length_srping')

end
