clear;
clc;

% 创建函数
F_external = 2;
delta_t = 1e-3; T = 1;
duty = 0.2; % 占空比
Time = 0:delta_t:T;
square_func = zeros(1, length(Time));
square_func(1, 1:round(length(Time)*duty)) = 1;
Force_ext = F_external * square_func;
plot(Time, Force_ext);