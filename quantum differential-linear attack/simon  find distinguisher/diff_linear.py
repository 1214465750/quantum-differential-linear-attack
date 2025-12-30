import gurobipy as gp
from gurobipy import GRB, quicksum, max_, min_
import numpy as np
import math
import time


class SimonCryptoAnalysis:
    def __init__(self, word_size=16, nr_diff=5, nr_diffLin=4, nr_lin=4):
        """
        SIMON密码差分-线性分析模型 - 完整修复版
        修复重点：
        1. log函数定义域保护（abs_corr > 0）
        2. 约束命名唯一性
        3. 边界冲突检测
        4. 数值稳定性增强
        """
        # === 1. 初始化模型和参数 ===
        self.model = gp.Model("SIMON_DiffLinear_Cryptanalysis")
        self.word_size = word_size
        self.nr_diff = nr_diff
        self.nr_diffLin = nr_diffLin
        self.nr_lin = nr_lin
        self.constraint_counter = 0  # 约束计数器，确保唯一命名

        # === 2. 声明所有变量 ===
        self._declare_variables()

        # === 3. 构建约束系统 ===
        self._build_differential_constraints()
        self._build_diff_linear_constraints()
        self._build_linear_constraints()

        # === 4. 设置目标函数 ===
        self._set_objective()

        # === 5. 添加保护机制 ===
        self._add_safeguards()

    def _generate_unique_name(self, base_name):
        """生成唯一约束名称"""
        self.constraint_counter += 1
        return f"{base_name}_{self.constraint_counter}_{int(time.time())}"

    def _declare_variables(self):
        """声明所有类型变量 - 修复键结构问题"""
        # 2.1 差分部分变量
        self.left_diff = self.model.addVars(self.nr_diff + 1, self.word_size, vtype=GRB.BINARY, name="left_diff")
        self.right_diff = self.model.addVars(self.nr_diff + 1, self.word_size, vtype=GRB.BINARY, name="right_diff")
        self.varibits = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="varibits")
        self.doublebits = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="doublebits")
        self.Belt = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="Belt")
        self.Z = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="Z")
        self.dummy_and = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="dummy_and")
        self.CC = self.model.addVars(self.nr_diff, self.word_size, vtype=GRB.BINARY, name="CC")
        self.diff_probability = self.model.addVar(vtype=GRB.INTEGER, name="diff_probability")

        # 2.2 差分-线性部分变量
        self.diffLin_left = self.model.addVars(
            self.nr_diffLin + 1, self.word_size, lb=-1, ub=1, name="diffLin_left")
        self.diffLin_right = self.model.addVars(
            self.nr_diffLin + 1, self.word_size, lb=-1, ub=1, name="diffLin_right")
        self.diffLin_corr = self.model.addVar(lb=-1, ub=1, name="diffLin_corr")

        # 关键修复：添加log函数保护变量
        self.abs_corr = self.model.addVar(lb=0, ub=1, name="abs_corr")
        self.abs_corr_safe = self.model.addVar(lb=1e-8, ub=1, name="abs_corr_safe")  # 保护变量
        self.diffLinComplement_corr = self.model.addVar(lb=-10, ub=0, name="diffLinComplement_corr")

        # 2.3 线性部分变量 - 修复键结构问题
        self.lin_state = {}
        for r in range(self.nr_lin + 1):
            for pos in [0, 1]:
                for i in range(self.word_size):
                    key = (r, pos, i)  # 使用三维键结构
                    self.lin_state[key] = self.model.addVar(
                        vtype=GRB.BINARY,
                        name=f"lin_state_{r}_{pos}_{i}"
                    )
        self.lin_and1 = self.model.addVars(self.nr_lin, self.word_size, vtype=GRB.BINARY, name="lin_and1")
        self.lin_and2 = self.model.addVars(self.nr_lin, self.word_size, vtype=GRB.BINARY, name="lin_and2")
        self.p_lin = self.model.addVars(self.nr_lin, self.word_size, vtype=GRB.BINARY, name="p_lin")
        self.lin_corr = self.model.addVar(vtype=GRB.INTEGER, name="lin_corr")

    def _add_safeguards(self):
        """添加数值保护机制 - 关键修复点"""
        # 1. LOG函数保护：确保输入值>0
        self.model.addGenConstrMax(
            self.abs_corr_safe,
            [self.abs_corr],
            1e-8,  # 确保最小值1e-8
            name=self._generate_unique_name("log_safeguard")
        )

        # 2. 边界一致性检查
        for v in self.model.getVars():
            if v.LB > v.UB:
                # 根据变量名重置边界
                if 'binary' in v.VarName:
                    v.setAttr("LB", 0)
                    v.setAttr("UB", 1)
                elif 'diffLin' in v.VarName:
                    v.setAttr("LB", -1)
                    v.setAttr("UB", 1)

    def _LRot(self, var_arr, shift):
        """循环左移实现"""
        return [var_arr[(i + shift) % self.word_size] for i in range(self.word_size)]

    def _continuous_xor_bit(self, x, y):
        """连续异或操作 (值域[-1,1])"""
        return -1.0 * x * y

    def _continuous_and_bit(self, x, y):
        """连续与操作 (值域[-1,1])"""
        return 0.25 * (x + y - x * y - 1.0)

    def _build_differential_constraints(self):
        """构建差分部分约束系统 - 修复运算符问题"""
        for rnd in range(self.nr_diff):
            # 获取当前轮输入
            x0 = [self.left_diff[rnd, i] for i in range(self.word_size)]
            x1 = [self.right_diff[rnd, i] for i in range(self.word_size)]

            # 中间计算：位移操作
            x2 = self._LRot(x0, 8)
            x3 = self._LRot(x0, 1)

            for i in range(self.word_size):
                # AND约束 - 使用唯一命名
                self.model.addGenConstrAnd(
                    self.dummy_and[rnd, i],
                    [x0[i], x2[i]],
                    name=self._generate_unique_name(f"AND_{rnd}_{i}")
                )

                # XOR约束（使用数学等价式替代^运算符）
                self.model.addConstr(
                    self.Belt[rnd, i] == self.dummy_and[rnd, i] + x1[i] - 2 * self.dummy_and[rnd, i] * x1[i],
                    name=self._generate_unique_name(f"XOR1_{rnd}_{i}")
                )

                # 另一个XOR约束
                self.model.addConstr(
                    self.Z[rnd, i] == self.Belt[rnd, i] + x3[i] - 2 * self.Belt[rnd, i] * x3[i],
                    name=self._generate_unique_name(f"XOR2_{rnd}_{i}")
                )

                # OR约束
                self.model.addGenConstrOr(
                    self.varibits[rnd, i],
                    [x0[i], x2[i]],
                    name=self._generate_unique_name(f"OR_{rnd}_{i}")
                )

                # AND约束
                self.model.addGenConstrAnd(
                    self.doublebits[rnd, i],
                    [x0[i], x2[i]],
                    name=self._generate_unique_name(f"AND_{rnd}_{i}")
                )

                # 状态转移 - 使用唯一命名
                self.model.addConstr(
                    self.left_diff[rnd + 1, i] == self.right_diff[rnd, i],
                    name=self._generate_unique_name(f"state_transfer1_{rnd}_{i}")
                )
                self.model.addConstr(
                    self.right_diff[rnd, i] == self.Z[rnd, i],
                    name=self._generate_unique_name(f"state_transfer2_{rnd}_{i}")
                )

        # 差分概率总和
        self.model.addConstr(
            self.diff_probability == quicksum(
                self.varibits[r, i] + 2 * self.doublebits[r, i]
                for r in range(self.nr_diff)
                for i in range(self.word_size)
            ),
            name=self._generate_unique_name("diff_prob_sum")
        )

    def _build_diff_linear_constraints(self):
        """构建差分-线性部分约束系统 - 修复运算符问题"""
        # 连接差分输出到差分-线性输入
        for i in range(self.word_size):
            # 二进制转连续变量 (0->-1, 1->1)
            self.model.addConstr(
                self.diffLin_left[0, i] == 2 * self.left_diff[self.nr_diff, i] - 1,
                name=self._generate_unique_name(f"cast_left_{i}")
            )
            self.model.addConstr(
                self.diffLin_right[0, i] == 2 * self.right_diff[self.nr_diff, i] - 1,
                name=self._generate_unique_name(f"cast_right_{i}")
            )

        # 差分-线性传播约束
        for rnd in range(self.nr_diffLin):
            for i in range(self.word_size):
                # 获取当前状态
                left = [self.diffLin_left[rnd, j] for j in range(self.word_size)]
                right = [self.diffLin_right[rnd, j] for j in range(self.word_size)]

                # SIMON轮函数近似
                rot_left = self._LRot(left, 8)
                rot_left2 = self._LRot(left, 1)

                # 连续与操作
                and_result = self.model.addVars(self.word_size, lb=-1, ub=1, name=f"and_res_{rnd}_{i}")
                for j in range(self.word_size):
                    expr = self._continuous_and_bit(left[j], rot_left[j])
                    self.model.addConstr(
                        and_result[j] == expr,
                        name=self._generate_unique_name(f"cont_and_{rnd}_{j}")
                    )

                # 异或操作（使用数学等价式替代^运算符）
                xor_result = self.model.addVars(self.word_size, lb=-1, ub=1, name=f"xor_res_{rnd}_{i}")
                for j in range(self.word_size):
                    expr = self._continuous_xor_bit(and_result[j], rot_left2[j])
                    self.model.addConstr(
                        xor_result[j] == expr,
                        name=self._generate_unique_name(f"cont_xor_{rnd}_{j}")
                    )

                # 更新状态
                for j in range(self.word_size):
                    self.model.addConstr(
                        self.diffLin_left[rnd + 1, j] == right[j],
                        name=self._generate_unique_name(f"state_update1_{rnd}_{j}")
                    )
                    self.model.addConstr(
                        self.diffLin_right[rnd + 1, j] == -1 * xor_result[j] * right[j],
                        name=self._generate_unique_name(f"state_update2_{rnd}_{j}")
                    )

        # 相关性计算 (log2(|r|))
        self.model.addGenConstrAbs(
            self.abs_corr,
            self.diffLin_corr,
            name=self._generate_unique_name("abs_corr")
        )

        # 关键修复：使用保护变量替代原始abs_corr
        self.model.addGenConstrLog(
            self.diffLinComplement_corr,
            self.abs_corr_safe,  # 使用保护变量
            name=self._generate_unique_name("log_corr_safe")
        )

    def _build_linear_constraints(self):
        """构建线性部分约束系统 - 修复键结构问题"""
        # 连接差分-线性输出到线性输入
        for i in range(self.word_size):
            key = (0, 0, i)  # 三维键
            # 连续变量转二进制 ([-1,1] -> {0,1})
            self.model.addConstr(
                self.lin_state[key] == (self.diffLin_left[self.nr_diffLin, i] + 1) / 2,
                name=self._generate_unique_name(f"cast_to_bin_{i}")
            )

            key = (0, 1, i)
            self.model.addConstr(
                self.lin_state[key] == (self.diffLin_right[self.nr_diffLin, i] + 1) / 2,
                name=self._generate_unique_name(f"cast_to_bin_right_{i}")
            )

        # 线性轮约束
        for rnd in range(self.nr_lin):
            for i in range(self.word_size):
                # 获取当前状态键
                state_key = (rnd, 0, i)
                state_key_shift8 = (rnd, 0, (i + 8) % self.word_size)
                state_key_shift1 = (rnd, 0, (i + 1) % self.word_size)
                state_key_right = (rnd, 1, i)

                # AND操作分解
                # 第一层: and1 = state0[i] & state0[(i+8)%ws]
                self.model.addGenConstrAnd(
                    self.lin_and1[rnd, i],
                    [self.lin_state[state_key], self.lin_state[state_key_shift8]],
                    name=self._generate_unique_name(f"lin_and1_{rnd}_{i}")
                )

                # 第二层: and2 = and1 & state0[(i+1)%ws]
                self.model.addGenConstrAnd(
                    self.lin_and2[rnd, i],
                    [self.lin_and1[rnd, i], self.lin_state[state_key_shift1]],
                    name=self._generate_unique_name(f"lin_and2_{rnd}_{i}")
                )

                # 状态更新
                next_key_left = (rnd + 1, 0, i)
                self.model.addConstr(
                    self.lin_state[next_key_left] == self.lin_state[state_key_right],
                    name=self._generate_unique_name(f"state_update_left_{rnd}_{i}")
                )

                next_key_right = (rnd + 1, 1, i)
                # 修复：使用数学等价式替代^运算符
                self.model.addConstr(
                    self.lin_state[next_key_right] ==
                    self.lin_and2[rnd, i] + self.lin_state[state_key] - 2 * self.lin_and2[rnd, i] * self.lin_state[
                        state_key],
                    name=self._generate_unique_name(f"state_update_right_{rnd}_{i}")
                )

                # 线性相关性
                self.model.addConstr(
                    self.p_lin[rnd, i] == self.lin_and2[rnd, i],
                    name=self._generate_unique_name(f"lin_corr_bit_{rnd}_{i}")
                )

        # 线性相关性总和
        self.model.addConstr(
            self.lin_corr == quicksum(
                self.p_lin[r, i]
                for r in range(self.nr_lin)
                for i in range(self.word_size)
            ),
            name=self._generate_unique_name("lin_corr_sum")
        )

    def _set_objective(self):
        """设置复合目标函数"""
        # 最小化: 差分概率 + 线性相关性 - 对数相关性
        total_cost = (
                self.diff_probability +
                self.lin_corr -
                self.diffLinComplement_corr
        )
        self.model.setObjective(total_cost, GRB.MINIMIZE)

    def solve(self, time_limit=7200):
        """求解模型"""
        try:
            # 设置求解参数
            self.model.Params.TimeLimit = time_limit
            self.model.Params.MIPFocus = 2  # 聚焦最优边界
            self.model.Params.IntegralityFocus = 1  # 强化整数解
            self.model.Params.NonConvex = 2  # 允许非凸二次约束
            self.model.Params.LogToConsole = 1  # 启用求解日志
            self.model.Params.NumericFocus = 3  # 提高数值稳定性
            self.model.Params.FeasibilityTol = 1e-9  # 更严格的可行性容差

            # 阶段验证（调试用）
            if not self._validate_phases():
                return False

            # 导出模型用于调试
            self.model.write("simon_model.lp")

            # 优化求解
            self.model.optimize()

            # 结果分析
            if self.model.status == GRB.OPTIMAL:
                self._print_solution()
                return True
            else:
                print(f"求解状态: {self.model.status}")
                if self.model.status == GRB.INFEASIBLE:
                    print("模型不可行，进行冲突分析...")
                    self._analyze_infeasibility()
                return False
        except gp.GurobiError as e:
            print(f"Gurobi 错误 [{e.errno}]: {e.message}")
        except Exception as e:
            print(f"未知错误: {str(e)}")
        return False

    def _validate_phases(self):
        """分阶段验证模型可行性"""
        phases = [
            ("差分约束", self._build_differential_constraints),
            ("差分-线性约束", self._build_diff_linear_constraints),
            ("线性约束", self._build_linear_constraints)
        ]

        # 临时模型用于阶段验证
        temp_model = gp.Model("temp_validation")
        temp_model.Params.LogToConsole = 0

        for phase_name, phase_func in phases:
            try:
                # 重置并添加当前阶段约束
                temp_model.remove(temp_model.getConstrs())
                phase_func()  # 在当前模型添加约束
                temp_model.update()
                temp_model.optimize()
                if temp_model.status == GRB.INFEASIBLE:
                    print(f"! 阶段 '{phase_name}' 导致不可行")
                    return False
            except Exception as e:
                print(f"阶段验证错误 ({phase_name}): {str(e)}")
                return False
        return True

    def _analyze_infeasibility(self):
        """详细分析不可行原因"""
        self.model.computeIIS()  # 冲突约束分析
        self.model.write("conflict.ilp")

        # 提取冲突详细信息
        print("=== 冲突约束分析 ===")
        for c in self.model.getConstrs():
            if c.IISConstr:
                print(f"冲突约束: {c.ConstrName} (类型: {c.Sense}, RHS={c.RHS})")

        # 提取冲突边界
        for v in self.model.getVars():
            if v.IISLB: print(f"冲突下界: {v.VarName} ≥ {v.LB}")
            if v.IISUB: print(f"冲突上界: {v.VarName} ≤ {v.UB}")

        print("导出冲突详情: conflict.ilp")

        # 尝试可行性松弛
        print("尝试可行性松弛...")
        relax_obj = self.model.feasRelaxS(
            relaxobjtype=0,  # 最小化惩罚
            minrelax=False,
            constrpenn=[1.0] * len(self.model.getConstrs())
        )
        self.model.optimize()
        print(f"松弛成本: {relax_obj}")

        # 分析显著松弛的约束
        for i, constr in enumerate(self.model.getConstrs()):
            slack = constr.getAttr('QCSlack')
            if slack > 0.1:  # 显著松弛
                print(f"高冲突约束: {constr.ConstrName} (松弛量={slack})")

    def _print_solution(self):
        """输出最优解"""
        print(f"\n=== 最优目标值: {self.model.ObjVal:.4f} ===")

        # 输出差分特征
        print("\n差分路径:")
        for r in range(self.nr_diff + 1):
            left_val = ''.join(str(int(self.left_diff[r, i].X)) for i in range(self.word_size))
            right_val = ''.join(str(int(self.right_diff[r, i].X)) for i in range(self.word_size))
            print(f"轮次 {r}: L={left_val}, R={right_val}")

        # 输出相关性
        print(f"\n差分-线性相关性: {self.diffLin_corr.X:.6f}")
        print(f"对数相关性: {self.diffLinComplement_corr.X:.6f}")

        # 保存完整模型
        self.model.write("simon_solution.sol")

# ======================
# 执行分析
# ======================
if __name__ == "__main__":
    try:
        print("启动 SIMON 密码差分-线性分析...")
        analyzer = SimonCryptoAnalysis(
            word_size=16,
            nr_diff=5,
            nr_diffLin=4,
            nr_lin=4
        )
        print("模型构建完成，开始求解...")
        success = analyzer.solve(time_limit=7200)  # 2小时求解时限
        if success:
            print("求解成功!")
        else:
            print("求解未达到最优解")
    except Exception as e:
        print(f"程序异常: {str(e)}")

