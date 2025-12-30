import cupy as cp
import numpy as np
from math import *
from time import time
half_block=16
m=4

z=cp.zeros(62,dtype=cp.bool_)
z=[1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,]


def binary_ones_parity():
    # 步骤 1: 生成查找表（0~65535 的每个值的二进制中1的奇偶性）
    lookup_table = cp.zeros(65536, dtype=cp.uint8)  # 初始化查找表
    for i in range(65536):
        lookup_table[i] = bin(i).count('1') % 2  # 奇数个1为1，偶数个为0
    lookup_table=lookup_table.astype(cp.bool_)
    # 步骤 2: 将输入数组转换为 uint16 类型（确保值在0~65535范围内）

    # 步骤 3: 利用查找表直接映射结果
    return lookup_table

def generate_mk(NK):
    # 生成索引数组 [0, 1, 2, ..., N-1]

    rng = np.random.default_rng()
    m1 = rng.choice(2 ** 32, size=NK, replace=False, ).astype(np.uint32)
    m2 = rng.choice(2 ** 32, size=NK, replace=False, ).astype(np.uint32)
    mk = np.empty((NK, 4), dtype=np.uint16)
    k0 = ((m1 >> 16) & 0xFFFF).astype(np.uint16)
    k1 = (m1 & 0xFFFF).astype(np.uint16)
    k2 = ((m2 >> 16) & 0xFFFF).astype(np.uint16)
    k3 = (m2 & 0xFFFF).astype(np.uint16)
    mk[:, 0] = k0
    mk[:, 1] = k1
    mk[:, 2] = k2
    mk[:, 3] = k3
    mk=cp.array(mk).astype(cp.uint16)
    return mk


def generate_p(N):
    # 计算总块数（每块含2^16个元素）

    rng = np.random.default_rng()
    m1 = rng.choice(2 ** 32, size=N, replace=False, ).astype(np.uint32)
    p = np.empty((N, 2), dtype=np.uint16)
    p0 = ((m1 >> 16) & 0xFFFF).astype(np.uint16)
    p1 = (m1 & 0xFFFF).astype(np.uint16)
    p[:, 0] = p0
    p[:, 1] = p1
    p = cp.array(p).astype(cp.uint16)
    return p

def diff_enc(diff_l,diff_r,lamda_l, lamda_r,round,round_wrong,N,NK):
    # 奇偶性表
    lookup_table = binary_ones_parity()
    list_cor=[]
    #求ELP，测试大量密钥
    mk=generate_mk(NK)
    p=generate_p(N)

    for i in range(NK):
        rk = keysche(mk[i], round+round_wrong)
        rk_wrong=cp.random.randint(0,2**16,round_wrong,dtype=cp.uint16)
        #生成明文


        p1_l = p[:, 0]
        p1_r = p[:, 1]
        p2_l = p1_l ^ diff_l
        p2_r = p1_r ^ diff_r

        c1_l, c1_r = enc(p1_l, p1_r, rk, round+round_wrong)
        c1_l, c1_r = dec(c1_l, c1_r,rk_wrong,round_wrong)
        c2_l, c2_r = enc(p2_l, p2_r, rk, round+round_wrong)
        c2_l, c2_r = dec(c2_l, c2_r, rk_wrong, round_wrong)

        a=(lamda_l&(c1_l^c2_l))^(lamda_r&(c1_r^c2_r))
        a=lookup_table[a]

        cor=1-2*cp.sum(a)/N
        #转换数据类型GPU-CPU
        cor=cor.get()
        list_cor.append(cor)


    return list_cor

def calELP(list_cor):
    len_cor=len(list_cor)
    he=0
    for i in list_cor:
        he+=i*i
    he=he/len_cor
    return he

def enc(p_l,p_r,rk,round):

    for i in range(round):
        p_l,p_r=oneenc(p_l,p_r,rk[i])
    return p_l,p_r

def dec(p_l,p_r,rk,round):

    for i in range(round):
        p_l,p_r=onedec(p_l,p_r,rk[i])
    return p_l,p_r
def oneenc(p_l,p_r,rk):

    return rk^p_r^sl(p_l,2)^(sl(p_l,1)&sl(p_l,8)),p_l

def onedec(p_l,p_r,rk):
    return p_r,rk^p_l^sl(p_r,2)^(sl(p_r,1)&sl(p_r,8))
def sl(a,b):
    b=b%half_block

    return ((a << b) ^ (a >> (half_block - b))) & 0xffff

def keysche(mk,round):
    rk = cp.zeros(round, dtype=cp.uint16)
    rk[:m] = mk

    for i in range(m, round):
        tmp = sl(rk[i - 1], -3)
        if m == 4:
            tmp ^= rk[i - 3]
        tmp = tmp ^ sl(tmp, -1)
        rk[i] = ~rk[i - m] ^ tmp ^ z[(i - m) % 62] ^ 3
    return rk

if __name__=='__main__':
    round = 8
    diff_l = 0x10
    diff_r = 0x60
    lamda_l = 0x200
    lamda_r = 0x80
    round_wrong=4
    N=2**20
    NK=2**15
    t1=time()

    list_cor = diff_enc(diff_l, diff_r, lamda_l, lamda_r, round,round_wrong,N,NK)

    elp=calELP(list_cor)
    list_cor.append(elp)
    list_cor.append(N)
    list_cor.append(NK)
    # print(list_cor)
    np.save('simon_wrong_0x%x,0x%x_0x%x,0x%x_R%d_RW%d_N%d_NK%d_DKP.npy'%(diff_l,diff_r,lamda_l,lamda_r,round,round_wrong ,log2(N),log2(NK)), list_cor)


    t2 = time()
    print('时间',t2-t1)




