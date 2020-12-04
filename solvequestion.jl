E = 210e9;
A = 0.005;
args = PoleArgs(E,A); # 杆件参数
nodes = [
    (0.0,   0.0),
    (5.0,   7.0),
    (5.0,   0.0),
    (10.0,  7.0),
    (10.0,  0.0),
    (15.0,  0.0),
]; # 节点

poles = [
    (1,2),
    (1,3),
    (2,3),
    (2,4),
    (2,5),
    (3,5),
    (4,5),
    (4,6),
    (5,6),
]; # 杆件

truss = System(); # 创建系统
node(truss, nodes); # 添加节点 
pole(truss, poles, args); # 添加杆件
constraint(truss, 1, true, true); # 添加约束
constraint(truss, 6, true, true); # 添加约束
external_force(truss, 2, (20e3, 0.00)); # 添加外力
truss |> generate_K |> generate_P |> generate_C |> solve_δ |> caculate_ε; # 生成, 求解

δ_mm = 1000* truss.δ;
σ = E*truss.ε;