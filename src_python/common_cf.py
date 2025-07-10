# common_cf.py
import numpy as np
from scipy import constants

# 物理常数
c = constants.c  # 光速 (m/s)
G = constants.G  # 万有引力常数 (N m^2 / kg^2)
M_sun = 1.989e30  # 太阳质量 (kg)
Mpc = 1e6 * constants.parsec # 百万秒差距 (m)

# 宇宙学参数 (Planck 2018)
H0 = 67.4  # 哈勃常数 (km/s/Mpc)
h = H0 / 100
Omega_m = 0.315
Omega_lambda = 0.685

def distance(z):
    """
    根据红移z计算光度距离 (单位: m)
    使用数值积分计算。
    """
    from scipy.integrate import quad
    
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)
        
    # 从0到z积分
    integral, _ = quad(integrand, 0, z)
    
    # 光度距离 d_L = (1+z) * comoving_distance
    # comoving_distance = (c/H0_in_SI) * integral
    H0_in_SI = (H0 * 1000) / Mpc  # H0 in s^-1
    
    d_L = (1 + z) * (c / H0_in_SI) * integral
    return d_L

if __name__ == '__main__':
    # 测试
    z_test = 0.1
    d_L_meters = distance(z_test)
    d_L_Mpc = d_L_meters / Mpc
    print(f"红移 z = {z_test}对应的光度距离是: {d_L_Mpc:.4f} Mpc")
