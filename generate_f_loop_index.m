function index = generate_f_loop_index(U_entire_name)

    U_entire_name_str = num2str(U_entire_name);
    row = double(U_entire_name_str(:,1))-48; % double(a): get the ASCII of a, then minus 48 to get the num of string a
    col = double(U_entire_name_str(:,3))-48;
    direction = double(U_entire_name_str(:,4))-48;
    point = double(U_entire_name_str(:,5))-48;

    index = [row,col,direction,point];

end