# plot_figure.py
import matplotlib.pyplot as plt
import corner
import numpy as np

def plot_corner(sampler, param_keys, truths, burn_in=500, thin=15):
    """
    使用corner库绘制角图
    
    sampler: emcee Sampler对象
    param_keys: 参数标签列表
    truths: 真实参数值
    burn_in: 要丢弃的老化步数
    thin: 采样间隔
    """
    flat_samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    
    fig = corner.corner(
        flat_samples, 
        labels=param_keys, 
        truths=truths,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True, 
        title_kwargs={"fontsize": 12}
    )
    plt.show()

def plot_fisher_ellipse(fisher_matrix, param_keys, truths):
    """
    (功能保留) 绘制Fisher矩阵预测的误差椭圆
    这是一个更复杂的绘图，这里只提供一个框架
    """
    print("绘制Fisher矩阵椭圆的功能可以进一步实现。")
    # 这需要从协方差矩阵 (Fisher矩阵的逆) 中提取2x2子矩阵，
    # 并计算特征值和特征向量来定义椭圆。
    pass
