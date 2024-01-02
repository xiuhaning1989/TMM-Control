function data = process_data()

    data_folder = 'data/metadata';
    files = dir(data_folder);
    data_names = {files.name};
    data_names(1:2) = [];
    data_names = erase(data_names, '.mat');
    data_names = data_names';
    L = numel(data_names);

    % data_names_mat stores all the data names like '10 2 4 2 1', means: 10N, (2,4) Y Point 1
    data_names_mat = zeros(numel(data_names), 5);
    for i = 1:L
        temp = data_names{i};
        temp = split(temp, '_');

        for j = 1:5
            data_names_mat(i, j) = str2double(temp{j});
        end

    end

    % change all the lattice name into individual numbers
    [U_entire_name, ~, ~, ~, ~] = store_u_name_unit_cell(5, 5);
    U_entire_name_str = num2str(U_entire_name);
    row = double(U_entire_name_str(:, 1)) - 48; % double(a): get the ASCII of a, then minus 48 to get the num of string a
    col = double(U_entire_name_str(:, 3)) - 48;
    direction = double(U_entire_name_str(:, 4)) - 48;
    point = double(U_entire_name_str(:, 5)) - 48;

    % final data struct
    data = struct('F', {},'F_index', {},'max_x', {},'max_x_index', {},'max_y', {},'max_y_index', {},'max_norm', {},'max_norm_index', {});

    % Data Processing: traverse the data according the data_names and extract the required data
    for i = 1:L
        data_name = strcat('data/metadata/', data_names{i});
        load(data_name);
        data(i).F = data_names_mat(i,1);
        data(i).F_index = data_names_mat(i,2:end);
        U_displacement = U(1:2:end, 123); % HERE: extract the 123 row of the data(0.122s)

        % Delete upper-left red lattice data and lower-right blue lattice data, and we choose to assign the minimum value of the data as a "delete"
        % upper-left red lattice:   6.1x[1-3] => 131,132,133,134,135,136
        % lower-right blue lattice: 1.5x3 => 9,10, 2.5x1 => 35,36, 2.6x2 => 161,162
        [min_all, ~] = min(abs(U_displacement));
        U_displacement(131:136) = min_all;
        U_displacement(9:10) = min_all;
        U_displacement(35:36) = min_all;
        U_displacement(161:162) = min_all;

        % Pocess the displacement in the X direction
        % Get the max values and their index
        % Put them into the data struct
        U_displacement_x = U_displacement(1:2:end, :);
        [~, max_x_index] = max(abs(U_displacement_x));
        max_x = U_displacement_x(max_x_index);
        max_x_index = [row(max_x_index*2),col(max_x_index*2),point(max_x_index*2)];
        data(i).max_x = max_x;
        data(i).max_x_index = max_x_index;

        % Pocess the displacement in the Y direction
        U_displacement_y = U_displacement(2:2:end, :);
        [~, max_y_index] = max(abs(U_displacement_y));
        max_y = U_displacement_y(max_y_index);
        max_y_index = [row(max_y_index*2),col(max_y_index*2),point(max_y_index*2)];
        data(i).max_y = max_y;
        data(i).max_y_index = max_y_index;
        
        % Pocess the displacement
        U_displacement_norm = sqrt(U_displacement_x .^ 2 + U_displacement_y .^ 2);
        [~, max_norm_index] = max(U_displacement_norm);
        max_norm = U_displacement_norm(max_norm_index);
        max_norm_index = [row(max_norm_index*2),col(max_norm_index*2),point(max_norm_index*2)];
        data(i).max_norm = max_norm;
        data(i).max_norm_index = max_norm_index;
    end
    disp('The data was processed successfully!');
    save('data.mat', 'data');
end