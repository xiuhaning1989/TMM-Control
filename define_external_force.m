function [force_index, force_func] = define_external_force(force_index_num, F_external, Time)
    % Define the function of the external force and the lattice coordinates to which the external force is applied.
    % At present, the function is only applicable to apply a single directional force, and it is a return to zero square wave.
    % index: lattice coordinates where the external force is applied
    % F_external: Magnitude of external force

    row = force_index_num(1);
    col = force_index_num(2);
    direction = force_index_num(3);
    point = force_index_num(4);
    
    force_index = [];

    for j = 1:1
        index = strcat(num2str(row), num2str(col), num2str(direction), num2str(point));

        f_U_name1 = str2double(index);

        if j < 10
            f_U_name1 = f_U_name1 / 1e3;
        elseif j < 100
            f_U_name1 = f_U_name1 / 1e4;
        else
            f_U_name1 = f_U_name1 / 1e5;
        end

        force_index = [force_index; f_U_name1];
    end

    duty = 0.1; % 占空比
    square_func = zeros(1, length(Time));
    square_func(1, 1:round(length(Time) * duty)) = 1;
    force_func = F_external * square_func;

    % plot(Time, Force_ext);
end
