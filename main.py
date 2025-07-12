from src.simulation.parameters import SimulationParameters
from src.models.physical_model import PhysicalModel
from src.solvers.numerical_solver import NumericalSolver
from src.visualization.output import Output
import os

def main():
    # 初始化参数
    params = SimulationParameters()
    print("initialize success")
    # 创建物理模型和求解器
    model = PhysicalModel(params)
    solver = NumericalSolver(params, model)
    
    # 进行模拟计算
    print("Starting simulation...")
    results = solver.solve()
    print("Simulation completed!")
    
    # 输出结果
    if params.plot_results:
        print("Plotting results...")
        Output.plot_results(results)
    
    # 保存结果
    if params.save_results:
        print("Saving results...")
        if not os.path.exists('results'):
            os.makedirs('results')
        Output.save_results(results, 'results/simulation_results.npz')
        print("Results saved to 'results/simulation_results.npz'")

if __name__ == "__main__":
    main()