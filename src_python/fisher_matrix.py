# fisher_matrix.py
import numpy as np
from waveform import Waveform
from detector import Detector

def calculate_fisher_matrix(params, param_keys, detector, f_gw):
    """
    使用数值微分计算Fisher矩阵
    
    params: 包含所有参数的字典
    param_keys: 要对其进行微分的参数键的列表
    detector: Detector对象
    f_gw: 频率数组
    """
    n_params = len(param_keys)
    fisher = np.zeros((n_params, n_params))
    
    # 获取天线方向图
    theta, phi, psi = params['theta'], params['phi'], params['psi']
    Fp, Fc = detector.antenna_pattern(theta, phi, psi)
    
    def get_h(p_dict):
        wf = Waveform(p_dict)
        hp, hc = wf.get_waveform(f_gw)
        return Fp * hp + Fc * hc

    # 数值微分
    h_derivs = []
    for i, key in enumerate(param_keys):
        p_temp = params.copy()
        
        # 使用中心差分来计算导数
        # 选择一个小的步长 dp
        dp = p_temp[key] * 1e-6 if p_temp[key] != 0 else 1e-8
        
        p_temp[key] = params[key] + dp
        h_plus = get_h(p_temp)
        
        p_temp[key] = params[key] - dp
        h_minus = get_h(p_temp)
        
        dh_dp = (h_plus - h_minus) / (2 * dp)
        h_derivs.append(dh_dp)
        
    # 计算内积 (h1|h2) = 4 * real(integral(h1* h2_conj / Sh_f df))
    df = np.diff(f_gw)
    df = np.append(df, df[-1]) # 保持维度一致
    
    Sh_f = detector.Sh

    for i in range(n_params):
        for j in range(i, n_params):
            integrand = h_derivs[i] * np.conj(h_derivs[j]) / Sh_f
            integral = np.sum(4 * np.real(integrand) * df)
            fisher[i, j] = integral
            fisher[j, i] = integral # Fisher矩阵是对称的
            
    return fisher
