function data_save(data_external_force,data_max_displacement,U,U_filename)
    save(strcat('data/external_force/',U_filename), 'data_external_force');
    save(strcat('data/max_displacement/',U_filename),'data_max_displacement');
    save(strcat('data/metadata/',U_filename),'U');
end

