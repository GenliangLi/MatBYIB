# mcmc.py
import numpy as np
import emcee
from waveform import Waveform
from detector import Detector

class MCMC:
    def __init__(self, data, params_true, param_keys, param_priors, detector, f_gw):
        """
        data: 模拟的探测器数据 (h_obs)
        params_true: 真实参数字典
        param_keys: MCMC要采样的参数列表
        param_priors: 参数的先验范围 {'key': [min, max]}
        detector: Detector对象
        f_gw: 频率数组
        """
        self.data = data
        self.param_keys_to_sample = [k for k, v in param_priors.items() if v[0] != v[1]]
        self.fixed_params = {k: v[0] for k, v in param_priors.items() if v[0] == v[1]}
        
        self.priors = {k: v for k, v in param_priors.items() if v[0] != v[1]}
        self.detector = detector
        self.f_gw = f_gw
        self.ndim = len(self.param_keys_to_sample)
        self.params_true_array = np.array([params_true[key] for key in self.param_keys_to_sample])
        
    def _log_prior(self, theta):
        """对数先验"""
        for i, key in enumerate(self.param_keys_to_sample):
            min_val, max_val = self.priors[key]
            if not (min_val <= theta[i] <= max_val):
                return -np.inf
        return 0.0

    def _log_likelihood(self, theta):
        """对数似然函数"""
        # 将采样参数和固定参数合并
        params_dict = dict(zip(self.param_keys_to_sample, theta))
        params_dict.update(self.fixed_params)
        
        # 构造波形
        wf = Waveform(params_dict)
        hp, hc = wf.get_waveform(self.f_gw)
        
        # 应用天线方向图
        Fp, Fc = self.detector.antenna_pattern(params_dict['theta'], params_dict['phi'], params_dict['psi'])
        h_model = Fp * hp + Fc * hc
        
        # 计算似然
        residual = h_model - self.data
        Sh_f = np.maximum(self.detector.Sh, 1e-50) # 避免除以零
        df = np.diff(self.f_gw)
        df = np.append(df, df[-1])
        
        inner_product = np.sum(4 * np.real(residual * np.conj(residual) / Sh_f) * df)
        
        return -0.5 * inner_product

    def _log_probability(self, theta):
        """对数后验概率"""
        lp = self._log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self._log_likelihood(theta)

    def run_mcmc(self, nwalkers=50, nsteps=5000, progress=True):
        """运行MCMC采样"""
        # **关键改动**: 从先验分布中初始化walkers，而不是在一个小球里
        initial_state = np.zeros((nwalkers, self.ndim))
        for i, key in enumerate(self.param_keys_to_sample):
            min_val, max_val = self.priors[key]
            initial_state[:, i] = np.random.uniform(min_val, max_val, nwalkers)
            
        sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self._log_probability)
        sampler.run_mcmc(initial_state, nsteps, progress=progress)
        
        return sampler
