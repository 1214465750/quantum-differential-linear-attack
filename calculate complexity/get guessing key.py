import numpy as np
from pysat.solvers import Cadical103
from pysat.formula import IDPool, CNF
from time import time
def tuice_forward(round,lamdal,lamdar,a,b,c):
    dl=np.zeros((round+1,16),dtype=np.uint8)
    dr=np.zeros((round+1,16),dtype=np.uint8)
    sl = np.zeros((round+1, 16), dtype=np.uint8)
    sr = np.zeros((round+1, 16), dtype=np.uint8)

    for i in range(16):
        dl[0][i]=(lamdal>>i)%2
        dr[0][i]=(lamdar>>i)%2



    for i in range(round):
        for j in range(16):
            if dr[i][j]==1:
                dl[i][j - a] = 1
                dl[i][j - b] = 1
                dl[i][j - c] = 1
                dl[i + 1][j] = 1
                sl[i][j - a] = 1
                sl[i][j - b] = 1
            if sr[i][j]==1:
                sl[i][j - a] = 1
                sl[i][j - b] = 1
                sl[i][j - c] = 1
                sl[i + 1][j] = 1
        for j in range(16):
            if dl[i][j]==1:
                dr[i+1][j]=1
            if sl[i][j]==1:
                sr[i+1][j]=1


    # for i in range(round+1):
    #     print(dl[i][::-1],dr[i][::-1])
    # print()
    # for i in range(round+1):
    #     print(sl[i][::-1],sr[i][::-1])

    print(np.sum(sr[:-1]))


def tuice_backward(round,lamdal,lamdar,a,b,c):
    dl=np.zeros((round+1,16),dtype=np.uint8)
    dr=np.zeros((round+1,16),dtype=np.uint8)
    sl = np.zeros((round+1, 16), dtype=np.uint8)
    sr = np.zeros((round+1, 16), dtype=np.uint8)

    for i in range(16):
        dl[-1][i]=(lamdal>>i)%2
        dr[-1][i]=(lamdar>>i)%2



    for i in range(round,0,-1):
        for j in range(16):
            if dl[i][j]==1:
                dr[i][j - a] = 1
                dr[i][j - b] = 1
                dr[i][j - c] = 1
                dr[i - 1][j] = 1
                sr[i][j - a] = 1
                sr[i][j - b] = 1
            if sl[i][j]==1:
                sr[i][j - a] = 1
                sr[i][j - b] = 1
                sr[i][j - c] = 1
                sr[i - 1][j] = 1
        for j in range(16):
            if dr[i][j] == 1:
                dl[i - 1][j] = 1
            if sr[i][j] == 1:
                sl[i - 1][j] = 1

    for i in range(round+1):
        print(dl[i][::-1],dr[i][::-1])
    print()
    for i in range(round+1):
        print(sl[i][::-1],sr[i][::-1])

    print(np.sum(sl[1:]))
    print('lin=',np.sum(dl[0]+dr[0]))

if __name__ == '__main__':
    round=1
    lamdal=0xa
    lamdar = 0x1
    # a, b, c = 8, 1, 2
    a,b,c=5,0,1
    # tuice_forward(round,lamdal,lamdar,a,b,c)

    round = 5
    lamdal = 0x2
    lamdar = 0x5
    # a, b, c = 8, 1, 2
    a, b, c = 5, 0, 1
    tuice_backward(round,lamdal,lamdar,a,b,c)
