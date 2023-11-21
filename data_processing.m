function data_processing(U,U_entire_name)
    U_displacement = U(1:2:end,:);
    U_displacement_x = U_displacement(1:2:end,:);
    U_displacement_y = U_displacement(2:2:end,:);
    U_displacement_norm = sqrt(U_displacement_x.^2+U_displacement_y.^2);
    [m, i] = max(U_displacement_norm);
    U_entire_name_str = num2str(U_entire_name);
    row = double(U_entire_name_str(:,1))-48;
    col = double(U_entire_name_str(:,3))-48;
    direction = double(U_entire_name_str(:,4))-48;
    point = double(U_entire_name_str(:,5))-48;
    index_data = [row, col, direction,point, U_entire_name];
    disp(i);
end

