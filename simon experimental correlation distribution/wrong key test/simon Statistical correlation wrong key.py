import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import norm

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
    list_cor = np.load('simon_wrong_0x10,0x60_0x200,0x80_R8_RW4_N20_NK15_DKP.npy').tolist()
    NK=list_cor[-1]
    N=list_cor[-2]
    elp=list_cor[-3]
    list_cor=list_cor[:-3]
    print('wt(elp)', -math.log2(elp))
    len_cor=len(list_cor)


    # 绘制直方图与KDE

    # 输入数据（假设已存入列表data）
    data=list_cor

    B=(2**32-N)/(2**32-1)


    # 设置正态分布参数（用户自定义）


    # 创建画布
    plt.figure(figsize=(12, 10))
    # 绘制直方图
    counts, bins, _ = plt.hist(data, bins=200, density=True,alpha=0.5, color='lightgreen', edgecolor='white',label='数据分布')

    # 生成并绘制正态曲线


    mu = 0
    sigma2 = B/N+2**(-32)

    x = np.linspace(np.min(data), np.max(data), 200)
    pdf = norm.pdf(x, mu, math.sqrt(sigma2))
    plt.plot(x, pdf, 'c--', linewidth=1, label=r'$\mu = 0, \sigma^2 = B/N+2^{-n}$')



    # # 标记均值线
    # plt.axvline(mu, color='red', linestyle=':', linewidth=2, label='均值')

    # 标注与美化
    plt.title('Simon32/64,8轮,错误密钥,'+r'$\left( \left( 0x0,0x8 \right)\to \left( 0x104,0x0 \right) \right)$ '+r'$,N=2^{%d},NK=2^{%d}$'%(math.log2(N),math.log2(NK)), fontsize=14)
    plt.xlabel('数值', fontsize=12)
    plt.ylabel('概率密度', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(linestyle='--', alpha=0.4)

    # out_dir1 = 'D:\学生\马成栋\我的论文\量子差分线性攻击\图片\\'
    # out_dir2 = 'Simon32_64_8轮_4轮_错误密钥_(0x0,0x8)_(0x104,0x0)_N=2^{%d}_NK=2^{%d}.png' % (math.log2(N), math.log2(NK))
    # plt.savefig(out_dir1 + out_dir2)

    plt.show()

