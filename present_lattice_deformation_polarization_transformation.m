function present_lattice_deformation_polarization_transformation(filename, U, Time, Coor_unit_cell_x, Coor_unit_cell_y, Alpha, Gamma, Theta, n, m, rotation_kappa, i_alpha, U_entire_name)

    Coor_unit_cell_x_name = zeros(n + 1, m, 9);
    Coor_unit_cell_y_name = zeros(n + 1, m, 9);

    Coor_unit_cell_x0 = Coor_unit_cell_x;
    Coor_unit_cell_y0 = Coor_unit_cell_y;
    % x, y coordinates for each unit cell from C, A, B, C, A', B' C, P, to Q

    % Choose the angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344)
    alpha = Alpha(i_alpha);
    gamma = Gamma(i_alpha);
    theta = Theta(i_alpha);
    % alpha=Alphanon(i_alpha);
    % gamma=Gammanon(i_alpha);
    % theta=Thetanon(i_alpha);

    % Plot all unit cells to form a n by m lattice
    figure;
    x_min = [];
    x_max = [];
    y_min = []
    y_max = [];
    for i = 1:n + 1

        for j = 1:m

            for k = 1:9
                XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
            end

            XY_unit = rotation_kappa * XY_unit;
            if i == 1
                % hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
                hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 1);
                x_min = [x_min,min(XY_unit(1,4:7))];
                x_max = [x_max,max(XY_unit(1,4:7))];
                y_min = [y_min,min(XY_unit(2,4:7))];
                y_max = [y_max,max(XY_unit(2,4:7))];
            elseif i == n + 1
                % hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
                hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 1);
                x_min = [x_min,min(XY_unit(1,1:4))];
                x_max = [x_max,max(XY_unit(1,1:4))];
                y_min = [y_min,min(XY_unit(2,1:4))];
                y_max = [y_max,max(XY_unit(2,1:4))];
            else
                % hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
                hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 1);
                % hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
                hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 1);
                hold on; plot(XY_unit(1, 8:9), XY_unit(2, 8:9), 'g-', 'linewidth', 1.5)
                x_min = [x_min,min(XY_unit(1,:))];
                x_max = [x_max,max(XY_unit(1,:))];
                y_min = [y_min,min(XY_unit(2,:))];
                y_max = [y_max,max(XY_unit(2,:))];
                axis equal
            end

        end

    end
    % X Y limits
    x_lim = [min(x_min), max(x_max)];
    y_lim = [min(y_min), max(y_max)];

    % Round
    x_lim = [floor(min(x_min)), ceil(max(x_max))];
    y_lim = [floor(min(y_min)), ceil(max(y_max))];

    title(['Homogeneous Lattice of ' num2str(n) ' rows and ' ...
               num2str(m) ' columns with angle \alpha = ' num2str(alpha)])


    for i = 1:n + 1

        for j = 1:m
            % Store names of variables of a unit cell only including blue, on the bottom of thr lattice
            % u_ijx3, u_i+1jx1, u_i+1j+1x2, u_ijy3, u_i+1jy1, u_i+1j+1y2

            result = strcat(num2str(i), num2str(j), num2str(1), num2str(3));
            Coor_unit_cell_x_name(i, j, 1) = str2num(result);
            Coor_unit_cell_x_name(i, j, 4) = str2num(result);
            Coor_unit_cell_x_name(i, j, 7) = str2num(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(3));
            Coor_unit_cell_y_name(i, j, 1) = str2num(result);
            Coor_unit_cell_y_name(i, j, 4) = str2num(result);
            Coor_unit_cell_y_name(i, j, 7) = str2num(result);
            result = strcat(num2str(i), num2str(j), num2str(1), num2str(1));
            Coor_unit_cell_x_name(i, j, 2) = str2num(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(1));
            Coor_unit_cell_y_name(i, j, 2) = str2num(result);
            result = strcat(num2str(i), num2str(j), num2str(1), num2str(2));
            Coor_unit_cell_x_name(i, j, 3) = str2num(result);
            result = strcat(num2str(i), num2str(j), num2str(2), num2str(2));
            Coor_unit_cell_y_name(i, j, 3) = str2num(result);
            result = strcat(num2str(i + 1), num2str(j), num2str(1), num2str(1));
            Coor_unit_cell_x_name(i, j, 5) = str2num(result);
            result = strcat(num2str(i + 1), num2str(j), num2str(2), num2str(1));
            Coor_unit_cell_y_name(i, j, 5) = str2num(result);
            result = strcat(num2str(i + 1), num2str(j + 1), num2str(1), num2str(2));
            Coor_unit_cell_x_name(i, j, 6) = str2num(result);
            result = strcat(num2str(i + 1), num2str(j + 1), num2str(2), num2str(2));
            Coor_unit_cell_y_name(i, j, 6) = str2num(result);

            if j < 10
                Coor_unit_cell_x_name(i, j, :) = Coor_unit_cell_x_name(i, j, :) / 1e3;
                Coor_unit_cell_y_name(i, j, :) = Coor_unit_cell_y_name(i, j, :) / 1e3;
            elseif j < 100
                Coor_unit_cell_x_name(i, j, :) = Coor_unit_cell_x_name(i, j, :) / 1e4;
                Coor_unit_cell_y_name(i, j, :) = Coor_unit_cell_y_name(i, j, :) / 1e4;
            else
                Coor_unit_cell_x_name(i, j, :) = Coor_unit_cell_x_name(i, j, :) / 1e5;
                Coor_unit_cell_y_name(i, j, :) = Coor_unit_cell_y_name(i, j, :) / 1e5;
            end

        end

    end

    % gif
    interval = ceil(length(Time)/15);
    T = Time(1:interval:end);
    U = U(:, 1:end);
    x = 0; y = 0;
    figure;

    for tt = 1:length(T)
        time0 = find(abs(Time - T(tt)) < 1e-4);
        u_displacement = U(1:2:end, time0);

        for i = 1:n + 1

            for j = 1:m

                for k = 1:7
                    p_entire = find(U_entire_name == Coor_unit_cell_x_name(i, j, k));
                    q_entire = find(U_entire_name == Coor_unit_cell_y_name(i, j, k));

                    %find U_entire_name(q_entire)==Coor_unit_cell_x(y)_name(i,j,k)
                    if isempty(p_entire) == 0
                        Coor_unit_cell_x(i, j, k) = Coor_unit_cell_x0(i, j, k) + u_displacement(p_entire);
                    end

                    if isempty(q_entire) == 0
                        Coor_unit_cell_y(i, j, k) = Coor_unit_cell_y0(i, j, k) + u_displacement(q_entire);
                    end

                end

            end

        end

        shift_horizontal1 = [Coor_unit_cell_x(2, 2, 2) - Coor_unit_cell_x(2, 1, 2), ...
                                 Coor_unit_cell_y(2, 2, 2) - Coor_unit_cell_y(2, 1, 2)];
        angle1 = -atan2(shift_horizontal1(2), shift_horizontal1(1));
        rotation_kappa1 = [cos(angle1), -sin(angle1);
                           sin(angle1), cos(-angle1)];
        minX = []; minY = []; maxX = []; maxY = [];
        plot(x, y)

        for i = 1:n + 1

            for j = 1:m

                for k = 1:7
                    XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
                end

                XY_unit(1, 8) = (XY_unit(1, 1) + XY_unit(1, 2)) / 2; XY_unit(2, 8) = (XY_unit(2, 1) + XY_unit(2, 2)) / 2;
                XY_unit(1, 9) = (XY_unit(1, 4) + XY_unit(1, 5)) / 2; XY_unit(2, 9) = (XY_unit(2, 4) + XY_unit(2, 5)) / 2;
                % XY_unit=rotation_kappa*XY_unit;
                XY_unit = rotation_kappa1 * XY_unit;
                minX = [minX, min(XY_unit(1, 1:3))]; minY = [minY, min(XY_unit(2, 4:6))];
                maxX = [maxX, max(XY_unit(1, :))]; maxY = [maxY, max(XY_unit(2, :))];

                if i == 1
                    % hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
                    hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 1);
                elseif i == n + 1
                    % hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
                    hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 1);
                else
                    % hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
                    hold on; plot(XY_unit(1, 1:4), XY_unit(2, 1:4), 'r-', 'linewidth', 1);
                    % hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
                    hold on; plot(XY_unit(1, 4:7), XY_unit(2, 4:7), 'b-', 'linewidth', 1);
                    hold on; plot(XY_unit(1, 8:9), XY_unit(2, 8:9), 'k--', 'linewidth', 1.5)
                end

            end

        end

        axis equal
        hold off
        title(['Time = ' num2str(T(tt), '%.2f') ' s'])

        xlim(x_lim);
        ylim(y_lim);

        f1 = getframe(gcf);
        imind1 = frame2im(f1);
        [imind1, cm1] = rgb2ind(imind1, 256);

        if tt == 1
            imwrite(imind1, cm1, filename, 'gif', 'Loopcount', 500, 'DelayTime', 0);
        else
            imwrite(imind1, cm1, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
        end

    end

end
