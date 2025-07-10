# waveform.py
import numpy as np
from common_cf import *

class Waveform:
    def __init__(self, params):
        """
        params: 包含波形参数的字典
        - Mc: 啁啾质量 (kg)
        - q:  质量比 (m2/m1 <= 1)
        - dL: 光度距离 (m)
        - tc: 合并时间
        - phic: 合并相位
        - iota: 轨道倾角
        """
        self.Mc = params['Mc']
        self.q = params['q']
        self.dL = params['dL']
        self.tc = params['tc']
        self.phic = params['phic']
        self.iota = params['iota']
        
        # 从啁啾质量和质量比计算总质量和对称质量比
        self.eta = self.q / (1 + self.q)**2
        self.M_total = self.Mc / self.eta**(3/5)

    def get_waveform(self, f_gw):
        """
        计算给定频率下的引力波傅里叶域波形 h(f)
        使用stationary phase approximation (SPA)
        """
        # 计算PN系数
        v = (np.pi * self.M_total * G / c**3 * f_gw)**(1/3)
        
        # 振幅 A(f)
        A = (1 / self.dL) * np.sqrt(5 * np.pi / 24) * \
            (G * self.Mc / c**3)**(5/6) * \
            (np.pi * f_gw)**(-7/6)
            
        # 相位 Psi(f)
        # 这里我们只使用到2PN阶，和MATLAB代码保持一致
        psi_0 = 2 * np.pi * f_gw * self.tc - self.phic - np.pi/4
        psi_pn = (3 / (128 * self.eta * v**5)) * \
                 (1 + (20/9 + 20 * self.eta / 3) * v**2 - 16 * np.pi * v**3 + \
                  (10 * (3058673/1016064 + 5429 * self.eta / 504 + 617 * self.eta**2 / 72)) * v**4)
        
        Psi = psi_0 + psi_pn
        
        # 极化
        h_plus = A * (1 + np.cos(self.iota)**2) / 2 * np.exp(1j * Psi)
        h_cross = A * np.cos(self.iota) * np.exp(1j * (Psi + np.pi/2)) # hx is shifted by pi/2

        return h_plus, h_cross

if __name__ == '__main__':
    # 测试
    m1_msun = 1e6
    m2_msun = 5e5
    z = 0.5
    
    m1 = m1_msun * M_sun
    m2 = m2_msun * M_sun
    Mc_kg = ((m1 * m2)**(3/5)) / ((m1 + m2)**(1/5))
    q_val = m2 / m1
    dL_m = distance(z)

    test_params = {
        'Mc': Mc_kg,
        'q': q_val,
        'dL': dL_m,
        'tc': 0,
        'phic': 0,
        'iota': np.pi / 3
    }
    
    f_gw_test = np.logspace(-4, 0, 1000)
    
    wf = Waveform(test_params)
    hp, hc = wf.get_waveform(f_gw_test)
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.loglog(f_gw_test, np.abs(hp), label='h_plus amplitude')
    plt.loglog(f_gw_test, np.abs(hc), label='h_cross amplitude', linestyle='--')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Strain Amplitude')
    plt.title('Frequency Domain Waveform')
    plt.legend()
    plt.grid(True)
    plt.show()
