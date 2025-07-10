# detector.py
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
# common_cf is not directly used here, but good practice to keep if dependencies change
from common_cf import *

class Detector:
    def __init__(self, detector_name, f_gw, t_obs_years=4.0):
        self.name = detector_name
        self.t_obs = t_obs_years * 365 * 24 * 3600  # 观测时间 (s)
        self.f_gw = f_gw # GW 频率 (Hz)
        self.Sh = self._load_sensitivity_curve()

    def _load_sensitivity_curve(self):
        """从文件中加载探测器灵敏度曲线"""
        try:
            # 使用pandas读取数据，更稳健
            data = pd.read_csv(f"{self.name}.txt", delim_whitespace=True, header=None, names=['freq', 'asd'])
            # ASD (Amplitude Spectral Density) to PSD (Power Spectral Density) -> Sh = ASD^2
            interp_func = interp1d(data['freq'], data['asd']**2, kind='linear', fill_value="extrapolate")
            return interp_func(self.f_gw)
        except FileNotFoundError:
            raise FileNotFoundError(f"未找到灵敏度曲线文件: {self.name}.txt")

    def antenna_pattern(self, theta, phi, psi):
        """计算天线方向图 F+ 和 Fx"""
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        cos_psi = np.cos(psi)
        sin_psi = np.sin(psi)

        # 假设探测器在黄道面上，并且臂与x和y轴成45度角
        # 这是LISA/Taiji的简化模型
        F_plus = 0.5 * (1 + cos_theta**2) * np.cos(2 * phi) * np.cos(2 * psi) - cos_theta * np.sin(2 * phi) * np.sin(2 * psi)
        F_cross = 0.5 * (1 + cos_theta**2) * np.cos(2 * phi) * np.sin(2 * psi) + cos_theta * np.sin(2 * phi) * np.cos(2 * psi)

        return F_plus, F_cross

if __name__ == '__main__':
    # 测试
    f_test = np.logspace(-4, 0, 100)
    taiji_detector = Detector('TAIJISEN', f_test)
    
    import matplotlib.pyplot as plt
    plt.loglog(f_test, np.sqrt(taiji_detector.Sh))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Strain ASD (1/sqrt(Hz))')
    plt.title('TAIJI Sensitivity Curve')
    plt.grid(True)
    plt.show()

    Fp, Fc = taiji_detector.antenna_pattern(np.pi/3, np.pi/4, np.pi/6)
    print(f"天线方向图示例: F+ = {Fp:.4f}, Fx = {Fc:.4f}")
