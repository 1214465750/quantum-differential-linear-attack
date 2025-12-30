import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import norm
import os


plt.rcParams['font.family'] = 'sans-serif'  # 使用通用无衬线字体族
plt.rcParams['font.sans-serif'] = [
    'Microsoft YaHei',   # 微软雅黑 (Windows)
    'SimHei',            # 黑体 (Windows)
    'WenQuanYi Micro Hei', # 文泉驿 (Linux)
    'Songti SC',         # 宋体 (macOS)
    'DejaVu Sans'        # 英文回退字体
]
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常



if __name__=='__main__':

    # 创建画布
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    ax = axs[0, 0]
    list_cor = np.load('simon_right_0x10,0x60_0x200,0x80_R8_N20_NK15_DKP.npy').tolist()
    NK = list_cor[-1]
    N = list_cor[-2]
    elp = list_cor[-3]
    list_cor = list_cor[:-3]
    print('wt(elp)', -math.log2(elp))
    len_cor = len(list_cor)
    data = list_cor
    B = (2 ** 32 - N) / (2 ** 32 - 1)
    # 绘制直方图 - 使用ax.hist代替plt.hist
    counts, bins, _ = ax.hist(data, bins=50, density=True, alpha=0.5, color='lightgreen', edgecolor='white',label='statistical distribution')
    # 生成并绘制正态曲线
    mu = sum(list_cor) / NK
    sigma2 = elp - mu ** 2
    x = np.linspace(np.min(data), np.max(data), 200)
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'm--', linewidth=3, label=r'$\mu =c, \sigma^2 = ELP-c^2$')  # 使用ax.plot

    mu = sum(list_cor) / NK
    sigma2 = B / N + elp - mu ** 2
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'b--', linewidth=3, label=r'$\mu =c, \sigma^2 = B/N+ELP-c^2$')  # 使用ax.plot

    # 设置标题、标签等 - 使用ax.set_title, ax.set_xlabel等
    name = 'correct guessing key,' + r'$\left( \left( 0x10,0x60 \right)\to \left( 0x200,0x80 \right) \right)$ '
    ax.set_title(name, fontsize=14)
    ax.set_xlabel('Experimental correlation', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(linestyle='--', alpha=0.4)

    ax = axs[0, 1]
    list_cor = np.load('simon_right_0x40,0x4190_0x10,0x4_R8_N20_NK15_DKP.npy').tolist()
    NK = list_cor[-1]
    N = list_cor[-2]
    elp = list_cor[-3]
    list_cor = list_cor[:-3]
    print('wt(elp)', -math.log2(elp))
    len_cor = len(list_cor)
    data = list_cor
    B = (2 ** 32 - N) / (2 ** 32 - 1)
    # 绘制直方图 - 使用ax.hist代替plt.hist
    counts, bins, _ = ax.hist(data, bins=50, density=True, alpha=0.5, color='lightgreen', edgecolor='white',
                              label='statistical distribution')
    # 生成并绘制正态曲线
    mu = sum(list_cor) / NK
    sigma2 = elp - mu ** 2
    x = np.linspace(np.min(data), np.max(data), 200)
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'm--', linewidth=3, label=r'$\mu =c, \sigma^2 = ELP-c^2$')  # 使用ax.plot

    mu = sum(list_cor) / NK
    sigma2 = B / N + elp - mu ** 2
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'b--', linewidth=3, label=r'$\mu =c, \sigma^2 = B/N+ELP-c^2$')  # 使用ax.plot

    # 设置标题、标签等 - 使用ax.set_title, ax.set_xlabel等
    name = 'correct guessing key,' + r'$\left( \left( 0x40,0x4190 \right)\to \left( 0x10,0x4 \right) \right)$ '
    ax.set_title(name, fontsize=14)
    ax.set_xlabel('Experimental correlation', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(linestyle='--', alpha=0.4)

    ax = axs[1, 0]
    list_cor = np.load('simon_wrong_0x10,0x60_0x200,0x80_R8_RW4_N20_NK15_DKP.npy').tolist()
    NK = list_cor[-1]
    N = list_cor[-2]
    elp = list_cor[-3]
    list_cor = list_cor[:-3]
    print('wt(elp)', -math.log2(elp))
    len_cor = len(list_cor)
    data = list_cor
    q_low, q_high = np.quantile(data, [0.01, 0.99])  # 截取99%数据范围
    x = np.linspace(q_low, q_high, 200)
    B = (2 ** 32 - N) / (2 ** 32 - 1)
    # 绘制直方图 - 使用ax.hist代替plt.hist
    counts, bins, _ = ax.hist(data, bins=50, density=True, alpha=0.5, color='lightgreen', edgecolor='white',
                              label='statistical distribution',range=(q_low, q_high,))

    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    # 生成并绘制正态曲线
    mu = 0
    sigma2 = B / N + 2 ** (-32)
    # x = np.linspace(np.min(data), np.max(data), 200)
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'c--', linewidth=3, label=r'$\mu = 0, \sigma^2 = B/N+2^{-n}$')

    # 设置标题、标签等 - 使用ax.set_title, ax.set_xlabel等
    name = 'incorrect guessing key,' + r'$\left( \left( 0x10,0x60 \right)\to \left( 0x200,0x80 \right) \right)$ '
    ax.set_title(name, fontsize=14)
    ax.set_xlabel('Experimental correlation', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(linestyle='--', alpha=0.4)

    ax = axs[1, 1]
    list_cor = np.load('simon_wrong_0x40,0x4190_0x10,0x4_R8_RW4_N20_NK15_DKP.npy').tolist()
    NK = list_cor[-1]
    N = list_cor[-2]
    elp = list_cor[-3]
    list_cor = list_cor[:-3]
    print('wt(elp)', -math.log2(elp))
    len_cor = len(list_cor)
    data = list_cor
    q_low, q_high = np.quantile(data, [0.01, 0.99])  # 截取99%数据范围
    x = np.linspace(q_low, q_high, 200)
    B = (2 ** 32 - N) / (2 ** 32 - 1)
    # 绘制直方图 - 使用ax.hist代替plt.hist
    counts, bins, _ = ax.hist(data, bins=50, density=True, alpha=0.5, color='lightgreen', edgecolor='white',
                              label='statistical distribution', range=(q_low, q_high,))

    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    # 生成并绘制正态曲线
    mu = 0
    sigma2 = B / N + 2 ** (-32)
    # x = np.linspace(np.min(data), np.max(data), 200)
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    ax.plot(x, pdf, 'c--', linewidth=3, label=r'$\mu = 0, \sigma^2 = B/N+2^{-n}$')

    # 设置标题、标签等 - 使用ax.set_title, ax.set_xlabel等
    name = 'incorrect guessing key,' + r'$\left( \left( 0x40,0x4190 \right)\to \left( 0x10,0x4 \right) \right)$ '
    ax.set_title(name, fontsize=14)
    ax.set_xlabel('Experimental correlation', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(linestyle='--', alpha=0.4)

    # 如果是在多子图环境中，可能还需要调整子图之间的间距，可以使用plt.subplots_adjust或fig.tight_layout()
    # fig.tight_layout()
    # fig.suptitle(r'Simon32/64差分-线性区分器正确和错误猜测密钥的实验相关性统计',fontsize=16,  fontweight='bold', y=0.95,  x=0.5,  ha='center')

    plt.show()


    out_dir2 = r'D:\学生\马成栋\我的论文\量子差分线性攻击\图片\Simon32_64差分-线性区分器正确和错误猜测密钥的实验相关性.png'
    fig.savefig(out_dir2)



