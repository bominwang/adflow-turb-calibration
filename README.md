# ADflow — 湍流模型系数校准版

基于 [MDO Lab ADflow](https://github.com/mdolab/adflow) 的修改版本，支持 **SA 和 SST 湍流模型闭合系数的运行时修改**，用于贝叶斯校准与不确定性量化研究。

## 修改内容

原版 ADflow 将 SA 和 SST 闭合系数声明为 Fortran `parameter`（编译期常量），运行时不可修改。本分支移除了 `parameter` 属性，并通过 f2py Python 接口暴露全部系数。

### 修改文件（5 个）

| 文件 | 修改内容 |
|------|----------|
| `src/modules/paramTurb.F90` | SA + SST 常数改为可变；新增 4 个 setter 子程序 |
| `src/f2py/adflow.pyf` | 新增 `paramturb` 模块：22 个变量 + 4 个子程序 |
| `src/initFlow/initializeFlow.F90` | 确保 `referenceState` 中不调用 setter（避免覆盖用户系数） |
| `src/turbulence/SST.F90` | f1 混合函数中硬编码 0.09 替换为 `rSSTBetas` |
| `src/turbulence/turbUtils.F90` | 涡粘计算中硬编码 0.09 替换为 `rSSTBetas` |

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

## 构建与部署

本仓库**不提供预编译二进制**。正确的使用方式是：在已有 ADflow 编译环境中，用 `patch_adflow_turb.py` 补丁原版源码后重新编译。

> **注意**: 仓库中的 `adflow/libadflow.so` 是构建环境的产物，通过 git 克隆可能因行尾转换而损坏。请务必在目标环境中从源码编译。

### 前置条件

- 安装 [Docker Desktop](https://www.docker.com/products/docker-desktop/)，确保已启动

### 第一步：拉取 MDO Lab 官方 Docker 镜像

在 Windows PowerShell / CMD / Git Bash 中执行（不是在 Docker 容器内）：

```bash
docker pull mdolab/public:u22-gcc-ompi-stable
```

镜像约 5 GB，首次拉取需要等待。完成后会看到 `Status: Downloaded newer image` 或 `Status: Image is up to date`。

### 第二步：启动容器

```bash
docker run -it --name adflow-turb mdolab/public:u22-gcc-ompi-stable bash
```

成功后终端提示符变为容器内的用户：

```
mdolabuser@xxxxxxxxxxxx:~$
```

> **容器管理**：
> - 退出容器：输入 `exit`
> - 重新进入已停止的容器：`docker start -i adflow-turb`
> - 删除容器：`docker rm adflow-turb`
> - 如果容器名已被占用：换一个名字（如 `adflow-turb2`），或先 `docker rm` 旧容器

以下所有步骤均在容器内执行。

### 第三步：克隆本仓库

```bash
git clone https://github.com/bominwang/adflow-turb-calibration.git /tmp/repo
```

预期输出：

```
Cloning into '/tmp/repo'...
...
Receiving objects: 100% ...
```

### 第四步：运行补丁脚本

```bash
python3 /tmp/repo/patch_adflow_turb.py
```

预期输出（5 个文件依次被补丁）：

```
======================================================================
Patching ADflow for runtime SA + SST coefficient modification
======================================================================

[1/5] Patching paramTurb.F90 (SA + SST mutable with initialisers)...
  Backup: .../paramTurb.F90.bak
  Patched: .../paramTurb.F90

[2/5] Patching adflow.pyf (f2py interface)...
  Backup: .../adflow.pyf.bak
  Patched: .../adflow.pyf

[3/5] Patching initializeFlow.F90 (no setter calls in referenceState)...
  .../initializeFlow.F90: no setter calls found, OK.

[4/5] Patching SST.F90 (f1 blending: 0.09 -> rSSTBetas)...
  Patched: .../SST.F90

[5/5] Patching turbUtils.F90 (eddy viscosity: 0.09 -> rSSTBetas)...
  Patched: .../turbUtils.F90

======================================================================
Patch complete!
======================================================================
```

补丁脚本默认补丁路径为 `/home/mdolabuser/repos/adflow/src`（MDO Lab Docker 镜像中的标准位置）。如需指定其他路径：

```bash
ADFLOW_SRC=/your/adflow/src python3 /tmp/repo/patch_adflow_turb.py
```

### 第五步：重新编译 ADflow

```bash
cd /home/mdolabuser/repos/adflow
make clean && make
```

编译过程中会出现一些 `Warning: Type mismatch` 警告，这是 ADflow 原版代码的 CGNS 模块问题，**不影响使用**。只要没有 `Error` 导致编译中断即可。

编译成功的标志是最后几行出现：

```
Testing if module libadflow can be imported...
Module libadflow was successfully imported
```

然后安装 Python 包：

```bash
pip install -e .
```

预期输出：`Successfully installed adflow-2.12.1`

### 第六步：验证安装

#### 快速验证（系数接口）

需要一个网格文件。在 Windows 上**另开一个终端**，将网格复制进容器：

```bash
docker cp /path/to/adflow-turb-calibration/examples/NACA0012/mesh/n0012.cgns adflow-turb:/tmp/mesh.cgns
```

回到容器内终端，运行验证脚本：

```python
python3 << 'EOF'
from adflow import ADFLOW
import os
os.makedirs("/tmp/output", exist_ok=True)

solver = ADFLOW(options={
    "gridFile": "/tmp/mesh.cgns",
    "outputDirectory": "/tmp/output",
    "equationType": "RANS",
    "turbulenceModel": "SA",
    "nCycles": 1,
})
pt = solver.adflow.paramturb

# 检查默认值
print("=== SA defaults ===")
print(f"cb1={pt.rsacb1:.4f} (expect 0.1355)")
print(f"cw1={pt.rsacw1:.4f} (expect 3.2391)")

# 修改系数并验证 cw1 自动重算
pt.setsaconstants(0.14, 0.65, 0.7, 0.40, 7.0, 0.25, 1.8, 1.2, 0.5)
expected_cw1 = 0.14 / 0.40**2 + (1 + 0.65) / 0.7
print(f"\n=== After setsaconstants ===")
print(f"cb1={pt.rsacb1:.4f} (expect 0.14)")
print(f"cw1={pt.rsacw1:.4f} (expect {expected_cw1:.4f})")
print(f"cw1 match: {abs(pt.rsacw1 - expected_cw1) < 1e-6}")

# SST 系数
pt.setsstconstants(0.41, 0.35, 0.085, 0.85, 0.5, 0.075, 1.0, 0.856, 0.0828)
print(f"\n=== After setsstconstants ===")
print(f"a1={pt.rssta1:.4f} (expect 0.35)")
print(f"betas={pt.rsstbetas:.4f} (expect 0.085)")
print(f"\n=== ALL PASSED ===")
EOF
```

**安装正确的预期输出**：

```
=== SA defaults ===
cb1=0.1355 (expect 0.1355)
cw1=3.2391 (expect 3.2391)

=== After setsaconstants ===
cb1=0.1400 (expect 0.14)
cw1=3.2321 (expect 3.2321)
cw1 match: True

=== After setsstconstants ===
a1=0.3500 (expect 0.35)
betas=0.0850 (expect 0.085)

=== ALL PASSED ===
```

如果看到 `ALL PASSED` 和 `cw1 match: True`，说明 SA 和 SST 系数接口全部正常。

#### 完整算例验证（NACA 0012 Cp 分布）

挂载本仓库到容器后，可运行完整的 NACA 0012 验证算例：

```bash
# 先停止当前容器
exit

# 重新启动并挂载仓库目录（Windows 路径示例）
docker run -it --name adflow-turb-full \
  -v /path/to/adflow-turb-calibration:/workspace/repo \
  mdolab/public:u22-gcc-ompi-stable bash

# 在容器内重复第三～五步（patch + compile + install）
python3 /workspace/repo/patch_adflow_turb.py
cd /home/mdolabuser/repos/adflow && make clean && make && pip install -e .

# 运行 SA 验证（4 组算例：1 默认 + 3 随机系数）
cd /workspace/repo/examples/NACA0012
mpirun --oversubscribe -np 2 python3 run_naca0012_sa_verify.py
```

SA 模型使用默认系数时的预期结果：

| 指标 | 预期值 |
|------|--------|
| Cl | ≈ 0.287 |
| Cd | ≈ 0.013 |
| 收敛 | ~60 次迭代达到 L2 残差 1e-15 |

如果 4 组算例全部收敛，且不同系数对应的 Cl/Cd 和 Cp 分布有差异，说明系数修改在实际求解中生效。

SST 模型验证同理：

```bash
mpirun --oversubscribe -np 2 python3 run_naca0012_sst_verify.py
```

### 超算部署

在无法直接使用 Docker 的 HPC 环境中，可构建 Singularity/Apptainer 镜像：

```bash
# 方法 1: 直接从 Docker 镜像转换
singularity pull adflow-turb.sif docker://mdolab/public:u22-gcc-ompi-stable

# 方法 2: 使用 Dockerfile 构建自定义镜像后转换
docker build -t adflow-turb:latest .
singularity pull adflow-turb.sif docker-daemon://adflow-turb:latest
```

Dockerfile 示例：

```dockerfile
FROM mdolab/public:u22-gcc-ompi-stable

COPY patch_adflow_turb.py /tmp/patch_adflow_turb.py
RUN python3 /tmp/patch_adflow_turb.py \
    && cd /home/mdolabuser/repos/adflow \
    && make clean && make \
    && pip install -e .
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
