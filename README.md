# ADflow — 湍流模型系数校准版

基于 [MDO Lab ADflow](https://github.com/mdolab/adflow) 的修改版本，支持 **SA 和 SST 湍流模型闭合系数的运行时修改**，用于贝叶斯校准与不确定性量化研究。

## 修改内容

原版 ADflow 将 SA 和 SST 闭合系数声明为 Fortran `parameter`（编译期常量），运行时不可修改。本分支移除了 `parameter` 属性，并通过 f2py Python 接口暴露全部系数。

### 修改文件（3 个）

| 文件 | 修改内容 |
|------|----------|
| `src/modules/paramTurb.F90` | SA + SST 常数改为可变；新增 4 个 setter 子程序 |
| `src/f2py/adflow.pyf` | 新增 `paramturb` 模块：22 个变量 + 4 个子程序 |
| `src/initFlow/initializeFlow.F90` | 新增 `setSADefaults()` 和 `setSSTDefaults()` 初始化调用 |

### 暴露的变量

**SA 模型** — 13 个变量（9 个校准参数 + 4 个辅助参数）：

| 变量名 | 参数 | 默认值 | 说明 |
|--------|------|--------|------|
| `rsacb1` | c_b1 | 0.1355 | 产生项系数 |
| `rsacb2` | c_b2 | 0.622 | 扩散项系数 |
| `rsacb3` | **sigma** | 0.6667 | **注意：存储的是 sigma (=2/3)，不是 c_b3** |
| `rsak` | kappa | 0.41 | von Karman 常数 |
| `rsacv1` | c_v1 | 7.1 | 壁面阻尼系数 |
| `rsacw1` | c_w1 | 3.2391 | 派生量：cb1/k² + (1+cb2)/sigma |
| `rsacw2` | c_w2 | 0.3 | 耗散项系数 |
| `rsacw3` | c_w3 | 2.0 | 耗散项系数 |
| `rsact1` | c_t1 | 1.0 | 转捩函数 |
| `rsact2` | c_t2 | 2.0 | 转捩函数 |
| `rsact3` | c_t3 | 1.2 | 转捩函数 |
| `rsact4` | c_t4 | 0.5 | 转捩函数 |
| `rsacrot` | c_rot | 2.0 | 旋转修正 |

**SST 模型** — 9 个独立变量：

| 变量名 | 参数 | 默认值 |
|--------|------|--------|
| `rsstk` | kappa | 0.41 |
| `rssta1` | a1 | 0.31 |
| `rsstbetas` | beta* | 0.09 |
| `rsstsigk1` | sigma_k1 | 0.85 |
| `rsstsigw1` | sigma_w1 | 0.5 |
| `rsstbeta1` | beta_1 | 0.075 |
| `rsstsigk2` | sigma_k2 | 1.0 |
| `rsstsigw2` | sigma_w2 | 0.856 |
| `rsstbeta2` | beta_2 | 0.0828 |

### 子程序

| 子程序 | 说明 |
|--------|------|
| `setsadefaults()` | 重置全部 13 个 SA 变量为默认值 |
| `setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)` | 设置 9 个 SA 校准参数，自动重算 c_w1 |
| `setsstdefaults()` | 重置全部 9 个 SST 变量为默认值 |
| `setsstconstants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)` | 设置 9 个 SST 参数 |

## 使用方法

```python
from adflow import ADFLOW

solver = ADFLOW(options=aeroOptions)
pt = solver.adflow.paramturb

# --- SA: 设置 9 个校准参数 (c_w1 自动重算) ---
pt.setsaconstants(
    0.1355,    # cb1
    0.622,     # cb2
    2.0/3.0,   # sigma (存储在 rsacb3 中)
    0.41,      # kappa
    7.1,       # cv1
    0.3,       # cw2
    2.0,       # cw3
    1.2,       # ct3
    0.5        # ct4
)

# --- SST: 设置 9 个独立参数 ---
pt.setsstconstants(
    0.41,      # sstk
    0.31,      # a1
    0.09,      # betas
    0.85,      # sigk1
    0.5,       # sigw1
    0.075,     # beta1
    1.0,       # sigk2
    0.856,     # sigw2
    0.0828     # beta2
)

# --- 也可单独设置 ---
pt.rsacw2 = 0.055       # SA
pt.rsstbeta1 = 0.06     # SST
```

## 预编译二进制

仓库中包含预编译的 `libadflow.so`（位于 `adflow/libadflow.so`），构建环境：
- Docker 镜像：`mdolab/public:u22-gcc-ompi-stable`
- Ubuntu 22.04, GCC, OpenMPI

使用预编译二进制时，运行环境必须与构建环境一致（与 MDO Lab Docker 镜像相同的操作系统、MPI 和 Python 版本）。

## 从源码构建

### 使用补丁脚本

如果已有原版 ADflow 源码：

```bash
python patch_adflow_turb.py   # 补丁 3 个源文件（幂等操作）
cd /path/to/adflow
make clean && make
pip install -e .
```

### 直接使用本仓库

```bash
git clone https://github.com/bominwang/adflow-turb-calibration.git
cd adflow-turb-calibration
make clean && make
pip install -e .
```

## 验证

### 单元测试

全部 7 项测试通过：

```
discover            : PASS   (22 variables + 4 subroutines)
sa_defaults         : PASS
sst_defaults        : PASS
sa_readwrite        : PASS   (13/13)
sst_readwrite       : PASS   (9/9)
sa_setter           : PASS   (9 params + c_w1 auto-recompute)
sst_setter          : PASS   (9 params)
```

### NACA 0012 算例验证

在 NACA 0012 翼型上进行了系数敏感性验证（M=0.75, α=1.5°, 海拔 10000 m）。每个湍流模型运行 4 组算例：1 组默认系数 + 3 组随机采样系数，提取表面压力系数 (Cp) 分布进行对比。

**SA 模型**（7 个校准参数独立采样，ct3/ct4 保持默认）：

| 算例 | Cl | Cd |
|------|--------|---------|
| default | 0.2867 | 0.01326 |
| random_1 | 0.2864 | 0.01324 |
| random_2 | 0.2891 | 0.01324 |
| random_3 | 0.2882 | 0.01253 |

![SA Cp 对比](doc/naca0012_sa_cp_compare.png)

**SST 模型**（9 个校准参数全部独立采样）：

| 算例 | Cl | Cd |
|------|--------|---------|
| default | 0.2942 | 0.01789 |
| random_1 | 0.2965 | 0.01667 |
| random_2 | 0.2960 | 0.01435 |
| random_3 | 0.2886 | 0.01727 |

![SST Cp 对比](doc/naca0012_sst_cp_compare.png)

算例脚本位于 `examples/NACA0012/` 目录下。

## 上游版本

基于 [mdolab/adflow](https://github.com/mdolab/adflow) v2.12.1。

## 许可证

与上游一致：GNU Lesser General Public License (LGPL), version 2.1。详见 [LICENSE.md](LICENSE.md)。
