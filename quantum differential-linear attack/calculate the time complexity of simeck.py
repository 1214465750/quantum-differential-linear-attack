import math


def calculate_T_Q2(B, N, ELP, K_union, n, p, k,kin,kout):

    base = 2 ** K_union  # 2^{|K_union|}
    term1 = B / N  # B/N

    # base * term1 + (base - 1) / (2 ** n) + ELP
    denominator_inner = (base-1) * (term1+2**(-n)) + ELP

    # term1 + ELP
    numerator_inner = term1 + ELP

    #  4 * sqrt(numerator_inner / denominator_inner)
    ratio = numerator_inner / denominator_inner
    denom = 2 * math.sqrt(ratio)

    #  π / denom
    T_Q2_front = math.pi / denom

    #  4 / sqrt(2^n) * (2^(2p+1) * sqrt(2^(2p+1))) / sqrt(base) * sqrt(denominator_inner)
    #  2^(2p+1) * sqrt(2^(2p+1)) = 2^(3p + 1.5)
    #  4 * 2^(3p + 1.5 - n/2 - K_union/2) * sqrt(denominator_inner)
    exponent_part1 = 3 * p + 1.5 - (n / 2) - (K_union / 2)
    denom_part1 = 2 * (2 ** exponent_part1) * math.sqrt(denominator_inner)
    part1 = math.pi / denom_part1

    # (π / 4) * sqrt(2^(k - K_union))
    part2 = (math.pi / 4) * math.sqrt(2 ** (min(k-kin+kout,k+kin-kout)))


    inner_sum = part1 + part2


    T_Q2 = T_Q2_front * inner_sum

    return T_Q2



if __name__ == "__main__":


    ELP = 2 ** (-23.28)
    lin = 35
    kin =31
    kout = 2
    K_union=kin+kout
    p=12.02
    c = 2 ** (-p)
    k = 64
    n=32
    N=2**(math.ceil(2*p)+1) if (2*p-lin+1 >=0) else 2**(p+0.5*(1+lin))
    N=2**31
    print(math.log2(N))
    B=(2**n-N)/(2**n-1)







    result = calculate_T_Q2(B, N, ELP, K_union, n, p, k,kin,kout)
    print(f"计算得到的 T_Q2 = {result,math.log2(result)}")
    result=2*N+result
    print(f"计算得到的 T_Q1 = {result, math.log2(result)}")
    print('穷举',32+math.log2(math.pi/4))