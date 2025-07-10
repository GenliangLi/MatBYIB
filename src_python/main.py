# main.py
import numpy as np
import time
from common_cf import *
from detector import Detector
from waveform import Waveform
from fisher_matrix import calculate_fisher_matrix
from mcmc import MCMC
from plot_figure import plot_corner

def load_input_params(filename='input.txt'):
    """
    从input.txt加载参数，按行顺序读取。
    这是为了匹配原始MATLAB脚本可能使用的fscanf行为。
    """
    param_keys_in_order = [
        'm1_msun', 'm2_msun', 'z', 'iota', 'psi',
        'theta', 'phi' 
        # 添加其他按顺序排列的参数键
    ]
    
    params = {}
    try:
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
            
            if len(lines) < len(param_keys_in_order):
                raise ValueError("input.txt 文件中的行数少于预期的参数数量。")

            for i, key in enumerate(param_keys_in_order):
                # 从行中移除逗号并转换为浮点数
                value_str = lines[i].replace(',', '')
                params[key] = float(value_str)
                
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{filename}'")
        raise
    except ValueError as e:
        print(f"错误: 处理输入文件时出错。 {e}")
        raise
        
    return params

def main():
    # 1. 加载输入参数
    try:
        inputs = load_input_params()
    except (FileNotFoundError, ValueError):
        print("由于输入文件错误，程序终止。")
        return

    # 2. 设置物理参数
    m1 = inputs['m1_msun'] * M_sun
    m2 = inputs['m2_msun'] * M_sun
    z = inputs['z']
    
    Mc = ((m1 * m2)**(3/5)) / ((m1 + m2)**(1/5))
    q = m2 / m1 if m1 > m2 else m1 / m2
    dL = distance(z)
    
    # 将 tc 和 phic 固定为0，简化问题
    tc_val, phic_val = 0.0, 0.0
    
    params_true = {
        'Mc': Mc, 'q': q, 'dL': dL, 'tc': tc_val, 'phic': phic_val,
        'iota': inputs['iota'], 'psi': inputs['psi'],
        'theta': inputs['theta'], 'phi': inputs['phi']
    }
    
    # 定义频率范围
    f_gw = np.logspace(-4, -1, 2000)
    
    # 3. 初始化探测器
    detector = Detector('TAIJISEN', f_gw)
    
    # 4. 生成"真实"的引力波信号 (作为模拟数据)
    wf_true = Waveform(params_true)
    hp_true, hc_true = wf_true.get_waveform(f_gw)
    
    Fp_true, Fc_true = detector.antenna_pattern(params_true['theta'], params_true['phi'], params_true['psi'])
    h_true = Fp_true * hp_true + Fc_true * hc_true
    
    # 添加高斯噪声
    df = np.diff(f_gw)
    df = np.append(df, df[-1])
    # 确保没有零或负的Sh值
    Sh_safe = np.maximum(detector.Sh, 1e-50) 
    sigma = np.sqrt(Sh_safe / (4 * df))
    noise = (np.random.normal(0, sigma) + 1j * np.random.normal(0, sigma)) / np.sqrt(2)
    
    data = h_true + noise

    # --- Fisher Matrix 分析 ---
    print("开始计算Fisher矩阵...")
    fisher_start_time = time.time()
    
    fisher_param_keys = ['Mc', 'q', 'iota', 'psi', 'theta', 'phi']
    fisher = calculate_fisher_matrix(params_true, fisher_param_keys, detector, f_gw)
    
    try:
        cov_matrix = np.linalg.inv(fisher)
        errors = np.sqrt(np.diag(cov_matrix))
        
        print(f"Fisher矩阵计算耗时: {time.time() - fisher_start_time:.2f} 秒")
        print("Fisher矩阵预测的1-sigma误差:")
        for i, key in enumerate(fisher_param_keys):
            print(f"  d({key}) = {errors[i]:.2e}")
            
    except np.linalg.LinAlgError:
        print("Fisher矩阵是奇异的，无法求逆。可能是参数间存在强关联或信噪比太低。")

    print("\n" + "="*30 + "\n")

    # --- MCMC 分析 ---
    print("开始运行MCMC...")
    mcmc_start_time = time.time()

    mcmc_param_keys = ['Mc', 'q', 'iota', 'psi', 'theta', 'phi']
    priors = {
        'Mc': [Mc * 0.9, Mc * 1.1],
        'q': [q * 0.8, 1.0], 
        'iota': [0, np.pi],
        'psi': [0, 2 * np.pi],
        'theta': [0, np.pi],
        'phi': [0, 2 * np.pi],
    }
    
    # 传递给MCMC的完整参数集
    mcmc_all_keys = list(params_true.keys())
    # 为MCMC的固定参数定义一个完整的先验字典
    full_priors = priors.copy()
    full_priors['dL'] = [dL, dL]
    full_priors['tc'] = [tc_val, tc_val]
    full_priors['phic'] = [phic_val, phic_val]

    mcmc_runner = MCMC(data, params_true, mcmc_all_keys, full_priors, detector, f_gw)
    
    sampler = mcmc_runner.run_mcmc(nwalkers=100, nsteps=4000, progress=True)
    
    print(f"MCMC运行耗时: {time.time() - mcmc_start_time:.2f} 秒")
    
    # 5. 绘制结果
    truths_for_plot = [params_true[key] for key in mcmc_param_keys]
    plot_corner(sampler, mcmc_param_keys, truths_for_plot, burn_in=500, thin=15)


if __name__ == '__main__':
    main()
