function max_data = data_processing(U,U_entire_name)
    % get the max displacement(norm) of each time and return as the following format
    % | row | col | point | max_displacement |
    % |  1  |  1  |   3   |    4.2343        |
    %  ...

    % get the norm of the displacement
    U_displacement = U(1:2:end,:);
    U_displacement_x = U_displacement(1:2:end,:);
    U_displacement_y = U_displacement(2:2:end,:);
    U_displacement_norm = sqrt(U_displacement_x.^2+U_displacement_y.^2);

    % get the max and max index of the displacement
    [m, i] = max(U_displacement_norm);

    % get the unit name
    U_entire_name_str = num2str(U_entire_name);
    row = double(U_entire_name_str(:,1))-48; % double(a): get the ASCII of a, then minus 48 to get the num of string a
    col = double(U_entire_name_str(:,3))-48;
    % direction = double(U_entire_name_str(:,4))-48;
    point = double(U_entire_name_str(:,5))-48;
    % displacement_data = [row, col, point, U_displacement_norm];
    
    max_data = [row(i*2),col(i*2),point(i*2),m.'];
end

