%% PreSettings
    close all;
    clear;
    clear -global;
    clc;

    addpath('functions');
    % rows and columns of unit cells.
    n = 5; m = 5;

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

    % Solve angles bistable lattice
    [Alpha, Gamma, Theta] = solve_angles_bistable_lattice(Lattice_config);

    % Present homogeneous hexagon
    [Coor_unit_cell_x, Coor_unit_cell_y, rotation_kappa] = present_homogeneous_hexagon(Lattice_config, Alpha, Gamma, Theta, n, m, i_alpha);

    % Store u name unit cell
    [U_entire_name, U_bottom_name, U_left_name, U_top_name, U_right_name] = store_u_name_unit_cell(n, m);

    % Generate the basic index_map
    index_map = generate_index_map(Coor_unit_cell_x, Coor_unit_cell_y, rotation_kappa, m, n);

    % Process the matadata for analyzing
    % the process_data() function takes about 30 seconds to finish
    file = 'data.mat';
    if exist(file, 'file') == 2
        disp('data.mat found!');
        load data.mat
    else
        disp('data.mat not found! Start to process the data');
        data = process_data();
    end


%% Drawing

    % Gradient colors setting
    % GRADIENT_COLORS: 6 gradient colors defined for drawing
    BASE_COLOR_1 = [20, 102, 127]/255;
    BASE_COLOR_2 = [226, 121, 35]/255;
    START_END_COLORS = [163, 200, 220, 036, 125, 169, 008, 072, 094;
                        241, 162, 169, 215, 053, 063, 165, 005, 021;
                        169, 214, 187, 069, 180, 128, 030, 142, 086;
                        211, 185, 216, 132, 082, 157, 096, 037, 130;
                        242, 231, 146, 237, 221, 049, 193, 174, 037;
                        246, 189, 148, 233, 111, 022, 188, 076, 000];
    GRADIENT_COLORS = zeros(29,3,6);
    for i = 1:6
        GRADIENT_COLORS(:,:,i) = [generate_gradient_colors(START_END_COLORS(i,1:3), START_END_COLORS(i,4:6), 28); START_END_COLORS(i,7:9)/255];
    end

    % draw F and its response
    % Each F has 29 different value
    F_index_to_draw = [3,1,1,1]; % the coordinates of the Force to be drawing (less than six)
    circle_size = linspace(100, 300, 29);
    circle_size = exp(linspace(-5, -2, 29))*1000;
    F_loop = [1:6, 8:2:20, 25:5:100];
    F_num = size(F_index_to_draw);

    
    % Figure 1: draw the max displacement
    figure('Position',[600, 50, 900, 900]);
    for i = 1:F_num
        draw_data = draw_response(F_index_to_draw(i,:), data, 1, index_map, F_loop, circle_size, GRADIENT_COLORS, i);
    end
    % Draw base cells 
    minZ = min([data.max_norm]);
    for i = 1:n + 1

        for j = 1:m

            for k = 1:9
                XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
            end

            XY_unit = rotation_kappa * XY_unit;

            if i == 1
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            elseif i == n + 1
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
            else
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            end

        end

    end
    xlim([-3, 7]); ylim([-1, 6]); % zlim([0, 6]);
    xlabel('X'); ylabel('Y'); zlabel('Max Displacement');
    grid on;
    view(-10, 60); % the view of the 3-d figure
    title('Max Displacement');
    hold off;

    % Figure 2: draw the max displacement of X direction
    figure('Position',[600, 50, 900, 900]);
    for i = 1:F_num
        draw_data = draw_response(F_index_to_draw(i,:), data, 2, index_map, F_loop, circle_size, GRADIENT_COLORS, i);
    end
    % Draw base cells 
    minZ = min([data.max_norm]);
    for i = 1:n + 1

        for j = 1:m

            for k = 1:9
                XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
            end

            XY_unit = rotation_kappa * XY_unit;

            if i == 1
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            elseif i == n + 1
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
            else
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            end

        end

    end
    xlim([-3, 7]); ylim([-1, 6]); % zlim([0, 6]);
    xlabel('X'); ylabel('Y'); zlabel('Max Displacement');
    grid on;
    view(-10, 60); % the view of the 3-d figure
    title('Max Displacement of X Direction');
    hold off;


    % Figure 3: draw the max displacement of Y direction
    figure('Position',[600, 50, 900, 900]);
    for i = 1:F_num
        draw_data = draw_response(F_index_to_draw(i,:), data, 3, index_map, F_loop, circle_size, GRADIENT_COLORS, i);
    end

    % Draw base cells 
    minZ = min([data.max_norm]);
    for i = 1:n + 1

        for j = 1:m

            for k = 1:9
                XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
            end

            XY_unit = rotation_kappa * XY_unit;

            if i == 1
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            elseif i == n + 1
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
            else
                plot3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, 'color', BASE_COLOR_2, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 1:4), XY_unit(2, 1:4), ones(1, 4)*minZ, BASE_COLOR_2, 'linestyle', 'none'); hold on;
                plot3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, 'color', BASE_COLOR_1, 'linewidth', 1); hold on;
                % fill3(XY_unit(1, 4:7), XY_unit(2, 4:7), ones(1, 4)*minZ, BASE_COLOR_1, 'linestyle', 'none'); hold on;
            end

        end

    end
    xlim([-3, 7]); ylim([-1, 6]); % zlim([0, 6]);
    xlabel('X'); ylabel('Y'); zlabel('Max Displacement');
    grid on;
    view(-10, 60); % the view of the 3-d figure
    title('Max Displacement of Y Direction');
    hold off;










    % Function: Draw the response of different forces
    function draw_data = draw_response(F_index, data, displacement_type, index_map, F_loop, circle_size, GRADIENT_COLORS, order)
        % displacement: 
        % 1 => draw max displacement
        % 2 => draw max displacement of x direction
        % 3 => draw max displacement of y direction

        if(~exist('order','var'))
            order = 2;  % the order of the gradient color
        end

        % draw response
        draw_data = zeros(29, 4); % store the max displacements and their index

        if displacement_type == 1 % draw the max displacement
            for F = 1:29
                coor_x = []; coor_y = [];
                for i = 1:numel(data)
                    if isequal(data(i).F_index, F_index) && isequal(data(i).F, F_loop(F))
                        for j = 1:numel(index_map)
                            if isequal(index_map(j).index, data(i).max_norm_index)
                                coor_x = index_map(j).x;
                                coor_y = index_map(j).y;
                                break;
                            end
                        end
                        break;
                    end
                end
                draw_data(F,:) = [F_loop(F), coor_x, coor_y, data(i).max_norm];
            end
        elseif displacement_type == 2 % draw the max displacement of X direction
            for F = 1:29
                coor_x = []; coor_y = [];
                for i = 1:numel(data)
                    if isequal(data(i).F_index, F_index) && isequal(data(i).F, F_loop(F))
                        for j = 1:numel(index_map)
                            if isequal(index_map(j).index, data(i).max_x_index)
                                coor_x = index_map(j).x;
                                coor_y = index_map(j).y;
                                break;
                            end
                        end
                        break;
                    end
                end
                draw_data(F,:) = [F_loop(F), coor_x, coor_y, data(i).max_x];
            end
        elseif displacement_type == 3 % draw the max displacement of Y direction
            for F = 1:29
                coor_x = []; coor_y = [];
                for i = 1:numel(data)
                    if isequal(data(i).F_index, F_index) && isequal(data(i).F, F_loop(F))
                        for j = 1:numel(index_map)
                            if isequal(index_map(j).index, data(i).max_y_index)
                                coor_x = index_map(j).x;
                                coor_y = index_map(j).y;
                                break;
                            end
                        end
                        break;
                    end
                end
                draw_data(F,:) = [F_loop(F), coor_x, coor_y, data(i).max_y];
            end
        else
            error('ERROR: Input the correct displacment-type: 1, 2, 3');
        end

        scatter3(draw_data(:,2)', draw_data(:,3)', draw_data(:,4)',circle_size, GRADIENT_COLORS(:,:,order)); hold on;
        scatter3(draw_data(:,2)', draw_data(:,3)', draw_data(:,4)',circle_size, GRADIENT_COLORS(:,:,order), "filled");  hold on;
        % line(draw_data(:,2)', draw_data(:,3)', draw_data(:,4)', 'Color', GRADIENT_COLORS(end,:,order), 'LineWidth', 1); hold on;

        % draw the position of F
        F_index_new = F_index;
        F_index_new(3) = []; 
        for j = 1:numel(index_map)
            if isequal(index_map(j).index, F_index_new)
                coor_x = index_map(j).x;
                coor_y = index_map(j).y;
                break;
            end
        end
        scatter3(coor_x, coor_y, min([data.max_norm]),100, GRADIENT_COLORS(end,:,order), "LineWidth",2);
        % Draw arrow
        if F_index(3) == 1
            quiver3(coor_x,coor_y,min([data.max_norm]),0.5,0,0,'Color', GRADIENT_COLORS(end,:,order), 'LineWidth', 2);
        elseif F_index(3) == 2
            quiver3(coor_x,coor_y,min([data.max_norm]),0,0.5,0,'Color', GRADIENT_COLORS(end,:,order), 'LineWidth', 2);
        end
        hold on; 
    end

    % Function: generate gradient colors according the starting and endding color
    function gradient_colors = generate_gradient_colors(start_color, end_color, num)
        step = (end_color - start_color) / (num - 1);
    
        gradient_colors = zeros(num, 3);
    
        for i = 1:num
            gradient_colors(i, :) = start_color + (i - 1) * step;
        end
    
        gradient_colors = gradient_colors / 255;
    end
    

