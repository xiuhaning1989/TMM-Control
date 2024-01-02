function index_map = generate_index_map(Coor_unit_cell_x, Coor_unit_cell_y, rotation_kappa, m, n)

    % Generate the index map of all unit cells
    % index => x_coor , y_coor
    % [1,1,3] => 0.3028, -0.3979

    temp = [];

    for i = 1:n + 1

            for j = 1:m

                for k = 1:9
                    XY_unit(1, k) = Coor_unit_cell_x(i, j, k); XY_unit(2, k) = Coor_unit_cell_y(i, j, k);
                end
                
                XY_unit = rotation_kappa * XY_unit;
                if i == 1
                    temp = [temp; i, j, 3, XY_unit(1,1), XY_unit(2,1)];
                else
                    temp = [temp; i, j, 1, XY_unit(1,2), XY_unit(2,2)];
                    temp = [temp; i, j, 2, XY_unit(1,3), XY_unit(2,3)];
                    temp = [temp; i, j, 3, XY_unit(1,1), XY_unit(2,1)];
                end

                if j == m
                    temp = [temp; i+1, j+1, 2, XY_unit(1,6), XY_unit(2,6)]; % right edge
                end

            end

    end

    temp(end,:) = [];

    index_map = struct('index', {}, 'x', {}, 'y', {});
    for i = 1:length(temp)
        index_map(i).index = temp(i,1:3);
        index_map(i).x = temp(i,4);
        index_map(i).y = temp(i,5);
    end

end