## 2023/11/16 修改记录

## 修改代码结构

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

## 一些具体调整

+ Lattice configuration改为结构体变量储存，以增加可读性

```matlab
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

```matlab
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


# 2023-11-23 修改记录

新增数据生成脚本，采用并行计算加速

## 新增数据生成脚本`Data_Generate_Main.m`

新增数据生成脚本`Data_Generate_Main.m`，根据`TMM_Main.m`修改而来。

### 具体说明

+ 循环最终采用两层，外层遍历所有的外力作用位置与方向（储存在`F_loop_index`变量中），内层遍历外力大小。

+ 采用MATLAB并行计算功能，需安装并行计算工具箱（Parallel Computing Toolbox），采用`parfor`对主循环进行并行。具体方法如下：

  ```matlab
  % 规定并行工作池数量（不大于CPU核心数）并开启并行计算
  parpool('local',6);
  
  % parfor代替for
  % 1. parfor无法嵌套使用
  % 2. 因为是将循环分割为多部分同时进行，故parfor循环中不能有迭代操作
  % 3. parfor循环变量i只能是整数
  parfor i = 1:N 
  	...
  end
  
  % 关闭并行工作池
  delete(gcp('nocreate'));
  ```

+ 理论上，并行`parfor`速度最大提升n倍，n为cpu核心数

### 新增函数

+ `generate_f_loop_index`: 用以生成矩阵形式储存的晶格位置名称用以循环
+ `data_processing`: 处理元数据`U`，得到所有时刻的最大位移以及对应位置坐标
+ `data_save`: 保存输出数据到`data`文件夹，由于`parfor`不支持在循环中直接调用`save`函数，故建立此函数用以保存数据。

## 修改文件结结构

新增functions文件夹存放所有函数, data文件夹存放生成的数据，具体目录结构图如下

```
├─data 存放生成数据
│  ├─external_force	存放外力大小、位置、方向
│  ├─max_displacement 存放最大位移
│  └─metadata 元数据
└─functions 函数
```







