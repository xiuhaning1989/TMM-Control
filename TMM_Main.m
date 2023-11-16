close all;
clear;
clear -global;
clc;

%% Params
    % rows and columns of unit cells.
    n = 4; m = 4;

    % Angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344).
    i_alpha = 16;

    % Lattice configuration.
    % Blue Triangle
    a_b = 0.5; b_b = 0.7; c_b = 1;
    psi_ab = acos((a_b ^ 2 - b_b ^ 2 - c_b ^ 2) / (-2 * b_b * c_b));
    psi_bb = acos((b_b ^ 2 - a_b ^ 2 - c_b ^ 2) / (-2 * a_b * c_b));
    psi_cb = acos((c_b ^ 2 - b_b ^ 2 - a_b ^ 2) / (-2 * b_b * a_b));

    % Red Triangle
    a_r = 0.4; b_r = 0.8; c_r = 1;
    psi_ar = acos((a_r ^ 2 - b_r ^ 2 - c_r ^ 2) / (-2 * b_r * c_r));
    psi_br = acos((b_r ^ 2 - a_r ^ 2 - c_r ^ 2) / (-2 * a_r * c_r));
    psi_cr = acos((c_r ^ 2 - b_r ^ 2 - a_r ^ 2) / (-2 * b_r * a_r));

    l_s = (a_b + b_r) / 2 * 0.95;
    k_bond = 1e4; k_spring = 100;
    eta_bond = 20; eta_spring = 10;
    m_b = 0.1; m_r = 0.1;

    Lattice_config.a_b = a_b;  
    Lattice_config.b_b = b_b;
    Lattice_config.c_b = c_b;
    Lattice_config.psi_ab = psi_ab;
    Lattice_config.psi_bb = psi_bb;
    Lattice_config.psi_cb = psi_cb;
    Lattice_config.a_r = a_r;
    Lattice_config.b_r = b_r;
    Lattice_config.c_r = c_r;
    Lattice_config.psi_ar = psi_ar;
    Lattice_config.psi_br = psi_br;
    Lattice_config.psi_cr = psi_cr;
    Lattice_config.l_s = l_s;
    Lattice_config.k_bond = k_bond;
    Lattice_config.k_spring = k_spring;
    Lattice_config.eta_bond = eta_bond;
    Lattice_config.eta_spring = eta_spring;
    Lattice_config.m_b = m_b;
    Lattice_config.m_r = m_r;

%% Main

    % Solve angles bistable lattice
    [Alpha, Gamma, Theta] = solve_angles_bistable_lattice(Lattice_config);

    % Present homogeneous hexagon
    [Coor_unit_cell_x, Coor_unit_cell_y,rotation_kappa] = present_homogeneous_hexagon(Lattice_config,Alpha,Gamma,Theta,n,m,i_alpha);

    % Store u name unit cell
    [U_entire_name, U_bottom_name, U_left_name, U_top_name, U_right_name] = store_u_name_unit_cell(n, m);

    % Solve dynamic equations 
    alpha = Alpha(i_alpha);
    gamma = Gamma(i_alpha);
    theta = Theta(i_alpha);

    theta_matrix = theta * ones(n + 1, m);
    alpha_matrix = alpha * ones(n + 1, m);

    % Solve dynamic equations
    % Define the y-direction of the bottom edge are always zero. Lattice can
    % slide at the bottom edge1
    j = 1;
    result = strcat(num2str(4), num2str(3), num2str(1), num2str(1)); U_fix_name(1, :) = str2double(result);
    result = strcat(num2str(4), num2str(3), num2str(2), num2str(1)); U_fix_name(2, :) = str2double(result);
    result = strcat(num2str(4), num2str(3), num2str(1), num2str(3)); U_fix_name(3, :) = str2double(result);
    result = strcat(num2str(4), num2str(3), num2str(2), num2str(3)); U_fix_name(4, :) = str2double(result);
    result = strcat(num2str(4), num2str(4), num2str(1), num2str(2)); U_fix_name(5, :) = str2double(result);
    result = strcat(num2str(4), num2str(4), num2str(2), num2str(2)); U_fix_name(6, :) = str2double(result);

    if j < 10
        U_fix_name = U_fix_name / 1e3;
    elseif j < 100
        U_fix_name = U_fix_name / 1e4;
    else
        U_fix_name = U_fix_name / 1e5;
    end

    % U_fix_name=[U_bottom_name(1:end);U_fix_name];

    % Define external force at the designated point of the lattice
    f_U_name = [];

    for j = 1:1
        result = strcat(num2str(3), num2str(1), num2str(1), num2str(2)); f_U_name1 = str2double(result);

        if j < 10
            f_U_name1 = f_U_name1 / 1e3;
        elseif j < 100
            f_U_name1 = f_U_name1 / 1e4;
        else
            f_U_name1 = f_U_name1 / 1e5;
        end

        f_U_name = [f_U_name; f_U_name1];
    end

    F_external = [1] .* ones(length(f_U_name), 1);

    delta_t = 1e-3; T = 1;
    Time = 0:delta_t:T;
    omega = 10 * pi;
    Force_ext = -F_external .* exp(-20 * (Time - 0.1) .^ 2) .* cos(omega * Time);

    % u=zeros(2*(6*n*m+2*n+2*m),length(Time));
    initialvals = zeros(2 * (6 * n * m + 2 * n + 2 * m), 1);

    tic;
    options2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-8 .* ones(1, 2 * (6 * n * m + 2 * n + 2 * m)));
    [t, U] = ode45(@(t, u) equation_motion_unit_cell_dissipation(t, u, n, m, ...
        alpha_matrix, theta_matrix, Lattice_config, U_entire_name, U_fix_name, ...
        Time, Force_ext, f_U_name), 0:delta_t:T, initialvals, options2);
    U = U';

    elapsedTime = toc;
    fprintf('Elapsed time: %.4f minutes\n', elapsedTime/60);
    
    % draw figure
    present_lattice_deformation_polarization_transformation('test1.gif', U, Time, Coor_unit_cell_x, Coor_unit_cell_y, Alpha, Gamma, Theta, n, m, rotation_kappa, i_alpha, U_entire_name)

