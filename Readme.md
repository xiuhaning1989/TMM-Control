## 2023/11/16 修改记录

### 修改代码结构

将所有Step进行封装，所有全局变量改为参数传递，代码结构如下：

`TMM_Main.m`

​	程序运行主函数，参数设置、每个Step运行都在这里。

`solve_angles_bistable_lattice.m`

​	原Step1 Solve angles bistable lattice的函数封装。

`present_homogeneous_hexagon.m`

​	原Step2 Present homogeneous hexagon的函数封装。

`store_u_name_unit_cell.m`

​	原Step3 Store u name unit cell的函数封装。

Step4没有封装，而是直接放在了主函数`TMM_Main.m`里面。

### 一些具体调整

+ Lattice configuration改为结构体变量储存，以增加可读性

```
Lattice_config.a_b = a_b;  
Lattice_config.b_b = b_b;
Lattice_config.c_b = c_b;
Lattice_config.psi_ab = psi_ab;
Lattice_config.psi_bb = psi_bb;
Lattice_config.psi_cb = psi_cb;
Lattice_config.a_r = a_r;
Lattice_config.b_r = b_r;
Lattice_config.c_r = c_r;
Lattice_config.psi_ar = psi_ar;
Lattice_config.psi_br = psi_br;
Lattice_config.psi_cr = psi_cr;
Lattice_config.l_s = l_s;
Lattice_config.k_bond = k_bond;
Lattice_config.k_spring = k_spring;
Lattice_config.eta_bond = eta_bond;
Lattice_config.eta_spring = eta_spring;
Lattice_config.m_b = m_b;
Lattice_config.m_r = m_r;
```

+ `str2num`函数替换为`str2double`，以增加些许性能。
+ 规范了一些代码格式。
+ 最后画图函数，添加了一个寻找XY坐标范围的功能，以固定绘制gif图时的坐标轴

```
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

```

最后`x_lim`和`y_lim`中储存了这次画图lattice的坐标上下限（已两头取整）