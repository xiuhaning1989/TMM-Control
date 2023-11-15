function [Coor_unit_cell_x, Coor_unit_cell_y,rotation_kappa] = present_homogeneous_hexagon(configuration, Alpha, Gamma, Theta, n, m, i_alpha)

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

    %Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
    % Coor_hexagon_x = ones(n - 1, m - 1, 7);
    % Coor_hexagon_y = ones(n - 1, m - 1, 7);

    % x, y coordinates for each hexagon from A, B, C, D, E, F, to A
    Coor_unit_cell_x = ones(n + 1, m, 9);
    Coor_unit_cell_y = ones(n + 1, m, 9);
    % x, y coordinates for each unit cell from C, A, B, C, A', B' C, P, to Q

    %Choose the angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344)
    alpha = Alpha(i_alpha);
    gamma = Gamma(i_alpha);
    theta = Theta(i_alpha);

    %%
    %Define transformation vectors of coordinates

    %From i,j to i+1,j+1
    shift_diagonal = [c_b - c_r * cos(gamma + psi_ab + psi_br) + a_b * cos(gamma - alpha + psi_ab + psi_br), ...
                        -c_r * sin(gamma + psi_ab + psi_br) + a_b * sin(gamma - alpha + psi_ab + psi_br)];

    %From i,1, to i+1,1
    kappa_vertical = 3 * pi - alpha - theta - psi_ar - psi_cr - psi_bb;
    shift_vertical = [b_r * cos(theta + psi_bb + psi_cr) + a_b * cos(kappa_vertical), ...
                        -b_r * sin(theta + psi_bb + psi_cr) + a_b * sin(kappa_vertical)];
    % phi=gamma+theta+psi_ab+psi_bb-2*pi; %phi=0 for homogeneous lattice

    %From 1,j to 1,j+1
    kappa_horizontal = psi_ab + gamma - pi;
    shift_horizontal = [c_b + a_r * cos(kappa_horizontal), a_r * sin(kappa_horizontal)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Represent an unit cell in a local coordinate system
    A = [0, 0];
    B = [-a_r * cos(psi_ab + gamma - pi), -a_r * sin(psi_ab + gamma - pi)];
    C = [b_r * cos(theta + psi_bb + psi_cr), -b_r * sin(theta + psi_bb + psi_cr)];
    Aa = [b_r * cos(theta + psi_bb + psi_cr) + a_b * cos(3 * pi - alpha - theta - psi_ar - psi_cr - psi_bb), ...
              -b_r * sin(theta + psi_bb + psi_cr) + a_b * sin(3 * pi - alpha - theta - psi_ar - psi_cr - psi_bb)];
    Bb = [b_r * cos(theta + psi_bb + psi_cr) - b_b * cos(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar), ...
              -b_r * sin(theta + psi_bb + psi_cr) + b_b * sin(theta + psi_bb + psi_cr + alpha + psi_cb + psi_ar)];
    P = [b_r / 2 * cos(theta + psi_bb + psi_cr), -b_r / 2 * sin(theta + psi_bb + psi_cr)];
    Q = [b_r * cos(theta + psi_bb + psi_cr) - a_b / 2 * cos(alpha + theta + psi_ar + psi_cr + psi_bb), ...
            -b_r * sin(theta + psi_bb + psi_cr) + a_b / 2 * sin(alpha + theta + psi_ar + psi_cr + psi_bb)];

    Coor_base_unit_cell = [C; A; B; C; Aa; Bb; C; P; Q] - shift_vertical;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reproduce n (row) by m (column) unit cells

    for j = 1:m % generate the 1st row of lattice (bottom edge)
        Coor_unit_cell_x(1, j, :) = Coor_base_unit_cell(:, 1) + (j - 1) * shift_horizontal(1);
        Coor_unit_cell_y(1, j, :) = Coor_base_unit_cell(:, 2) + (j - 1) * shift_horizontal(2);
    end

    for i = 1:n + 1 %generate the 1st column of the lattice (left edge)
        Coor_unit_cell_x(i, 1, :) = Coor_base_unit_cell(:, 1) + (i - 1) * shift_vertical(1);
        Coor_unit_cell_y(i, 1, :) = Coor_base_unit_cell(:, 2) + (i - 1) * shift_vertical(2);
    end

    for j = 1:m %generate unit cells diagonally from the bottom edge

        if m < n + 1

            for k = 0:m - j
                Coor_unit_cell_x(1 + k, j + k, :) = Coor_unit_cell_x(1, j, :) + k * shift_diagonal(1);
                Coor_unit_cell_y(1 + k, j + k, :) = Coor_unit_cell_y(1, j, :) + k * shift_diagonal(2);
            end

        else

            if j <= max(n + 1, m) - min(n + 1, m)

                for k = 0:min(n + 1, m) - 1
                    Coor_unit_cell_x(1 + k, j + k, :) = Coor_unit_cell_x(1, j, :) + k * shift_diagonal(1);
                    Coor_unit_cell_y(1 + k, j + k, :) = Coor_unit_cell_y(1, j, :) + k * shift_diagonal(2);
                end

            else

                for k = 0:max(n + 1, m) - j
                    Coor_unit_cell_x(1 + k, j + k, :) = Coor_unit_cell_x(1, j, :) + k * shift_diagonal(1);
                    Coor_unit_cell_y(1 + k, j + k, :) = Coor_unit_cell_y(1, j, :) + k * shift_diagonal(2);
                end

            end

        end

    end

    for i = 1:n + 1 %generate unit cells diagonally from the left edge

        if m < n + 1

            if i <= max(n + 1, m) - min(n + 1, m)

                for k = 0:min(n + 1, m) - 1
                    Coor_unit_cell_x(i + k, 1 + k, :) = Coor_unit_cell_x(i, 1, :) + k * shift_diagonal(1);
                    Coor_unit_cell_y(i + k, 1 + k, :) = Coor_unit_cell_y(i, 1, :) + k * shift_diagonal(2);
                end

            else

                for k = 0:max(n + 1, m) - i
                    Coor_unit_cell_x(i + k, 1 + k, :) = Coor_unit_cell_x(i, 1, :) + k * shift_diagonal(1);
                    Coor_unit_cell_y(i + k, 1 + k, :) = Coor_unit_cell_y(i, 1, :) + k * shift_diagonal(2);
                end

            end

        else

            for k = 0:n + 1 - i
                Coor_unit_cell_x(i + k, 1 + k, :) = Coor_unit_cell_x(i, 1, :) + k * shift_diagonal(1);
                Coor_unit_cell_y(i + k, 1 + k, :) = Coor_unit_cell_y(i, 1, :) + k * shift_diagonal(2);
            end

        end

    end

    %Plot all unit cells to form a n by m lattice
    angle = -atan2(shift_horizontal(2), shift_horizontal(1));
    rotation_kappa = [cos(angle), -sin(angle);
                    sin(angle), cos(-angle)];

    % figure;

    % for i = 1:n + 1

    %     for j = 1:m

    %         for k = 1:9
    %             XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
    %         end

    %         XY_unit = rotation_kappa * XY_unit;

    %         if i == 1
    %             hold on; fill(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b')
    %             hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 0.8);
    %         elseif i == n + 1
    %             hold on; fill(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r')
    %             hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 0.8);
    %         else
    %             hold on; fill(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r')
    %             hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 0.8);
    %             hold on; fill(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b')
    %             hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 0.8);
    %             hold on; plot(XY_unit(1, 8:9), XY_unit(2, 8:9), 'g-', 'linewidth', 1.5)
    %             axis equal
    %         end

    %     end

    % end

    % axis equal

    % title(['Homogeneous Lattice of ' num2str(n) ' rows and ' ...
    %            num2str(m) ' columns with angle \alpha = ' num2str(alpha)])
end
