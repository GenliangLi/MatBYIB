# common_cf.py
import numpy as np
from scipy import constants

# ==================================
# 1. 物理常数 (Physical Constants)
# ==================================
c = constants.c          # 光速 (m/s)
G = constants.G          # 万有引力常数 (N m^2 / kg^2)
M_sun = 1.989e30         # 太阳质量 (kg)
Mpc = 1e6 * constants.parsec # 百万秒差距 (m)

# =======================================
# 2. 宇宙学参数 (Cosmological Parameters)
# (Based on Planck 2018)
# =======================================
H0 = 67.4                 # 哈勃常数 (km/s/Mpc)
h = H0 / 100              # 无量纲哈勃常数

Omega_m = 0.315           # 物质密度参数
Omega_lambda = 0.685      # 暗能量密度参数

# 这个是被遗漏的常数
Omega_k = 1.0 - Omega_m - Omega_lambda  # 曲率密度参数 (在此模型中为 0)

# =======================================
# 3. 核心函数 (Core Function)
# =======================================
def distance(z):
    """
    根据红移z计算光度距离 (单位: m)
    使用数值积分计算。
    """
    from scipy.integrate import quad
    
    def integrand(z_prime):
        # 完整的积分项，包含了曲率项
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_k * (1 + z_prime)**2 + Omega_lambda)
        
    # 从0到z积分
    integral, _ = quad(integrand, 0, z)
    
    # 同移距离 (Comoving Distance)
    # H0_in_SI 是以 s^-1 为单位的哈勃常数
    H0_in_SI = (H0 * 1000) / Mpc
    d_C = (c / H0_in_SI) * integral
    
    # 光度距离 (Luminosity Distance)
    # d_L = (1+z) * d_C  (在平坦宇宙中)
    # 更通用的公式涉及到 d_M (横向同移距离)，但在平坦宇宙中 d_M = d_C
    d_L = (1 + z) * d_C
    
    return d_L

if __name__ == '__main__':
    # 测试
    print(f"定义的宇宙学参数:")
    print(f"  Omega_m = {Omega_m}")
    print(f"  Omega_lambda = {Omega_lambda}")
    print(f"  Omega_k = {Omega_k} (计算得出)")

    z_test = 0.1
    d_L_meters = distance(z_test)
    d_L_Mpc = d_L_meters / Mpc
    print(f"\n红移 z = {z_test}对应的光度距离是: {d_L_Mpc:.4f} Mpc")