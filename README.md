# Mg颗粒燃烧仿真模型

## 项目概述

本项目实现了Mg颗粒在CO2环境中的加热、熔融和反应过程的仿真模拟。模型采用分阶段求解策略，使用MATLAB内置的ODE求解器，避免复杂的自定义数值计算方法。

整个过程分为以下几个阶段：

1. **预热阶段**（T < 923K）：颗粒温度低于熔点，只有温度随时间变化，热量来源于辐射和对流换热。
2. **熔融阶段**（T = 923K）：颗粒达到熔点，温度保持不变，定义熔融百分数，所有热量用于内部金属Mg的熔化。
3. **加热阶段**（923K < T < 1366K）：熔融完成后，温度继续上升，直到达到沸点。
4. **气相燃烧阶段**（T ≥ 1366K）：颗粒达到沸点，发生气相燃烧（框架代码，留待后续实现）。

## 文件结构

```
Mg-CO2/
├── main.m                 # 主程序入口
├── Parameters.m           # 参数定义类
├── ParticleModel.m        # 颗粒物理模型类
├── HeatTransfer.m         # 热传递计算类
├── PreheatingStage.m      # 预热阶段求解器
├── MeltingStage.m         # 熔融阶段求解器
├── HeatingStage.m         # 熔融后加热阶段求解器
├── VaporizationStage.m    # 气相燃烧阶段求解器框架
├── StageManager.m         # 阶段管理器
├── Visualization.m        # 结果可视化
├── README.md              # 项目说明文档
└── results/               # 结果输出目录
    ├── simulation_results.mat  # MATLAB格式结果
    └── simulation_results.csv  # CSV格式结果
```

## 模块功能说明

### 1. 主程序 (`main.m`)

主程序入口，负责初始化参数、创建模型、运行仿真和可视化结果。

主要功能：
- 初始化参数
- 创建物理模型
- 创建阶段管理器
- 运行仿真
- 可视化和保存结果
- 输出仿真信息

### 2. 参数定义 (`Parameters.m`)

定义所有仿真参数的类，包括物理常数、材料属性、环境参数和仿真控制参数。

主要属性：
- 物理常数（气体常数、Stefan-Boltzmann常数等）
- 材料物理参数（Mg、MgO、C的密度、比热容等）
- 颗粒物理参数（初始直径、初始温度）
- 环境参数（环境温度、压力、发射率）
- 气体组成
- 仿真控制参数（时间步长、总时间）
- 求解器配置

主要方法：
- `Parameters()`: 构造函数
- `validate()`: 验证参数的有效性

### 3. 颗粒物理模型 (`ParticleModel.m`)

管理颗粒的物理状态和属性的类。

主要属性：
- 参数对象
- 热传递模型
- 颗粒状态（温度、直径、质量、熔融分数）

主要方法：
- `ParticleModel(params)`: 构造函数
- `calculateMass()`: 计算颗粒质量
- `determineStage()`: 确定当前阶段
- `calculateHeatCapacity()`: 计算颗粒的热容
- `calculateMeltingEnergy()`: 计算完全熔化所需的能量

### 4. 热传递计算 (`HeatTransfer.m`)

处理颗粒的热传递计算的类。

主要方法：
- `HeatTransfer(params)`: 构造函数
- `calculateConvection(T_particle, T_ambient, diameter)`: 计算对流换热
- `calculateRadiation(T_particle, T_ambient, diameter)`: 计算辐射换热
- `calculateTotalHeatTransfer(T_particle, T_ambient, diameter)`: 计算总热传递

### 5. 预热阶段求解器 (`PreheatingStage.m`)

处理颗粒从初始温度到熔点的加热过程。

主要方法：
- `PreheatingStage(model, params)`: 构造函数
- `solve(tspan)`: 求解预热阶段
- `odefun(t, y)`: ODE函数，定义温度变化率
- `events(t, y)`: 事件函数，当温度达到熔点时终止

### 6. 熔融阶段求解器 (`MeltingStage.m`)

处理颗粒在熔点温度下的熔融过程。

主要方法：
- `MeltingStage(model, params)`: 构造函数
- `solve(tspan)`: 求解熔融阶段
- `odefun(t, y)`: ODE函数，定义熔融分数变化率
- `events(t, y)`: 事件函数，当熔融分数达到1时终止

### 7. 熔融后加热阶段求解器 (`HeatingStage.m`)

处理颗粒从熔点到点火温度的加热过程。

主要方法：
- `HeatingStage(model, params)`: 构造函数
- `solve(tspan)`: 求解熔融后加热阶段
- `odefun(t, y)`: ODE函数，定义温度变化率
- `events(t, y)`: 事件函数，当温度达到点火温度时终止

### 8. 气相燃烧阶段求解器框架 (`VaporizationStage.m`)

处理颗粒在点火温度下的气相燃烧过程（框架代码，留待后续实现）。

主要方法：
- `VaporizationStage(model, params)`: 构造函数
- `solve(tspan)`: 求解气相燃烧阶段（框架）
- `odefun(t, y)`: ODE函数（待实现）
- `calculateReactionRate(T, ambient_composition)`: 计算反应速率（待实现）
- `calculateHeatRelease(reaction_rate)`: 计算反应放热（待实现）

### 9. 阶段管理器 (`StageManager.m`)

管理不同阶段的切换和求解。

主要属性：
- 颗粒模型
- 参数对象
- 各阶段求解器
- 结果存储（时间历史、温度历史、熔融分数历史、阶段历史）

主要方法：
- `StageManager(model, params)`: 构造函数
- `solve()`: 求解整个过程
- `recordResults(t, T, X_melt, stage_name)`: 记录结果

### 10. 结果可视化 (`Visualization.m`)

提供各种可视化方法。

主要方法：
- `plotResults(results)`: 绘制结果
- `saveResults(results, filename)`: 保存结果到文件
- `exportToCSV(results, filename)`: 导出结果到CSV文件

## 模拟流程

1. 初始化参数和模型
2. 确定初始阶段（通常是预热阶段）
3. 使用对应的阶段求解器求解
4. 检测阶段转换事件
5. 切换到下一个阶段
6. 重复步骤3-5，直到完成所有阶段或达到总仿真时间
7. 输出和可视化结果

## 物理模型

### 预热阶段

微分方程：
```
dT/dt = (Q_conv + Q_rad) / (m * Cp)
```
其中：
- Q_conv：对流换热
- Q_rad：辐射换热
- m：颗粒质量
- Cp：比热容

### 熔融阶段

微分方程：
```
dX_melt/dt = (Q_conv + Q_rad) / (m * L_fusion)
```
其中：
- X_melt：熔融分数（0-1）
- L_fusion：熔化潜热

### 熔融后加热阶段

微分方程：
```
dT/dt = (Q_conv + Q_rad) / (m * Cp)
```

### 气相燃烧阶段

待实现

## 使用方法

1. 在MATLAB环境中运行`main.m`
2. 结果将自动保存到`results`目录
3. 可视化结果将显示在图形窗口中

## 扩展方向

1. 完善气相燃烧阶段的实现
2. 添加更复杂的反应动力学模型
3. 考虑颗粒尺寸变化
4. 添加更多的环境参数影响
5. 实现三维可视化 