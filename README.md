# ADflow — 湍流模型系数校准版

基于 [MDO Lab ADflow](https://github.com/mdolab/adflow) 的修改版本，支持 **SA 和 SST 湍流模型闭合系数的运行时修改**，用于贝叶斯校准与不确定性量化研究。

## 原理

原版 ADflow 将 SA 和 SST 闭合系数声明为 Fortran `parameter`（编译期常量），运行时不可修改。本分支：

1. **移除 `parameter`**：将 `paramTurb.F90` 中的 SA 和 SST 系数改为可变模块变量，保留默认初始值
2. **添加 setter 子程序**：`setSAConstants`、`setSSTConstants` 等 Fortran 子程序
3. **替换硬编码常量**：`SST.F90` 和 `turbUtils.F90` 中的硬编码 `0.09`（beta*）替换为模块变量 `rSSTBetas`
4. **提供 ctypes Python API**：`adflow_turb_ctypes.py` 直接写入 `libadflow.so` 中 Fortran 符号的内存地址

## Python API

Python 接口通过 ctypes 直接访问 Fortran 模块符号，绕过 f2py（原因见[已知问题](#已知问题)）。

```python
from adflow import ADFLOW
from adflow_turb_ctypes import (
    set_sa_constants, set_sa_defaults, get_sa_constants,
    set_sst_constants, set_sst_defaults, get_sst_constants,
)

solver = ADFLOW(options=aeroOptions)

# --- SA: 设置 9 个校准参数 (cw1 自动重算) ---
set_sa_constants(
    cb1=0.1355, cb2=0.622, sigma=2.0/3.0, kappa=0.41,
    cv1=7.1, cw2=0.3, cw3=2.0, ct3=1.2, ct4=0.5,
)

# 读取当前值
print(get_sa_constants())
# {'cb1': 0.1355, 'cb2': 0.622, 'sigma': 0.6667, 'kappa': 0.41,
#  'cv1': 7.1, 'cw1': 3.2391, 'cw2': 0.3, 'cw3': 2.0, 'ct3': 1.2, 'ct4': 0.5}

# --- SST: 设置 9 个独立参数 ---
set_sst_constants(
    sstk=0.41, a1=0.31, betas=0.09,
    sigk1=0.85, sigw1=0.5, beta1=0.075,
    sigk2=1.0, sigw2=0.856, beta2=0.0828,
)

# 重置为默认值
set_sa_defaults()
set_sst_defaults()
```

### 使用注意事项

- **只创建一个求解器实例**，在不同 `AeroProblem` 之间切换系数。创建多个 `ADFLOW` 实例会导致 PETSc 崩溃。
- **每次求解前重新设置系数**：ADflow 内部切换 AeroProblem 时可能调用 `referenceState`。在 `solver(ap)` 之前设置系数。
- 将 `adflow_turb_ctypes.py` 复制到工作目录，或将本仓库根目录加入 `PYTHONPATH`。

### 系数参考

**SA 模型** — 9 个校准参数：

| 参数 | Fortran 变量 | 默认值 | 说明 |
|------|-------------|--------|------|
| cb1 | `rsaCb1` | 0.1355 | 产生项系数 |
| cb2 | `rsaCb2` | 0.622 | 扩散项系数 |
| sigma | `rsaCb3` | 2/3 | 扩散比（**注意：存储在 rsaCb3 中，不是 c_b3**） |
| kappa | `rsaK` | 0.41 | von Karman 常数 |
| cv1 | `rsaCv1` | 7.1 | 壁面阻尼系数 |
| cw2 | `rsaCw2` | 0.3 | 耗散项系数 |
| cw3 | `rsaCw3` | 2.0 | 耗散项系数 |
| ct3 | `rsaCt3` | 1.2 | 转捩函数系数 |
| ct4 | `rsaCt4` | 0.5 | 转捩函数系数 |

派生量：`cw1 = cb1/kappa^2 + (1+cb2)/sigma`（setter 自动重算）

辅助变量（可修改，但不在标准校准集中）：`rsaCt1` (1.0)、`rsaCt2` (2.0)、`rsaCrot` (2.0)

**SST 模型** — 9 个独立参数：

| 参数 | Fortran 变量 | 默认值 |
|------|-------------|--------|
| kappa | `rSSTK` | 0.41 |
| a1 | `rSSTA1` | 0.31 |
| beta* | `rSSTBetas` | 0.09 |
| sigma_k1 | `rSSTSigk1` | 0.85 |
| sigma_w1 | `rSSTSigw1` | 0.5 |
| beta_1 | `rSSTBeta1` | 0.075 |
| sigma_k2 | `rSSTSigk2` | 1.0 |
| sigma_w2 | `rSSTSigw2` | 0.856 |
| beta_2 | `rSSTBeta2` | 0.0828 |

## 构建与部署

本仓库**不提供预编译二进制**。工作流程：用 `patch_adflow_turb.py` 补丁 ADflow 源码，然后从源码编译。

> **注意**：仓库中的 `adflow/libadflow.so` 是构建环境的产物，通过 git 克隆可能因行尾转换而损坏。请务必在目标环境中从源码编译。

### 方式 A：Docker（推荐用于开发）

#### 1. 拉取 MDO Lab 官方 Docker 镜像

```bash
docker pull mdolab/public:u22-gcc-ompi-stable
```

#### 2. 启动容器

```bash
docker run -it --name adflow-turb mdolab/public:u22-gcc-ompi-stable bash
```

#### 3. 克隆并补丁

```bash
git clone https://github.com/bominwang/adflow-turb-calibration.git /tmp/repo
python3 /tmp/repo/patch_adflow_turb.py
```

预期输出（4 个文件被补丁）：

```
[1/4] Patching paramTurb.F90 (SA + SST mutable with initialisers)...
  Patched: .../paramTurb.F90
[2/4] Checking initializeFlow.F90 (no setter calls in referenceState)...
  .../initializeFlow.F90: no setter calls found, OK.
[3/4] Patching SST.F90 (f1 blending: 0.09 -> rSSTBetas)...
  Patched: .../SST.F90
[4/4] Patching turbUtils.F90 (eddy viscosity: 0.09 -> rSSTBetas)...
  Patched: .../turbUtils.F90
```

#### 4. 编译安装

```bash
cd /home/mdolabuser/repos/adflow
make clean && make
pip install -e .
```

编译过程中 `Type mismatch` 警告来自 ADflow 上游代码的 CGNS 模块，**不影响使用**。编译成功的标志：

```
Testing if module libadflow can be imported...
Module libadflow was successfully imported
```

#### 5. 验证

```bash
cp /tmp/repo/adflow_turb_ctypes.py /home/mdolabuser/
python3 -c "
import sys; sys.path.insert(0, '/home/mdolabuser')
from adflow_turb_ctypes import set_sa_constants, get_sa_constants, set_sa_defaults

set_sa_defaults()
c = get_sa_constants()
assert abs(c['cb1'] - 0.1355) < 1e-6, f'default failed: {c}'

set_sa_constants(0.5, 0.622, 2/3, 0.41, 7.1, 0.3, 2.0, 1.2, 0.5)
c = get_sa_constants()
assert abs(c['cb1'] - 0.5) < 1e-6, f'set failed: {c}'

set_sa_defaults()
c = get_sa_constants()
assert abs(c['cb1'] - 0.1355) < 1e-6, f'reset failed: {c}'

print('ctypes interface: OK')
print('SA read/write/reset: OK')
"
```

### 方式 B：HPC 原生编译

适用于不支持 Docker/Singularity 的超算环境。需要 GCC gfortran + OpenMPI。

可使用 `hpc/deploy.sh` 一键部署脚本，或按以下手动步骤操作。

#### 前置条件

- GCC gfortran（已测试 12.2）
- OpenMPI（已测试 4.1.5）或 Intel MPI
- PETSc 3.15+ （需编译为共享库）
- CGNS 4.x （需编译 Fortran 支持）
- Python 3 + numpy + mpi4py

#### 手动步骤

```bash
# 1. 克隆
git clone https://github.com/bominwang/adflow-turb-calibration.git
cd adflow-turb-calibration

# 2. 设置环境变量
export PETSC_DIR=/path/to/petsc-install
export PETSC_ARCH=""
export CGNS_HOME=/path/to/cgns-install

# 3. 补丁
export ADFLOW_SRC=$(pwd)/src
python3 patch_adflow_turb.py

# 4. 创建 config/config.mk（根据 HPC 环境调整路径）
cat > config/config.mk << 'EOF'
FF90 = mpifort
CC   = mpicc

FF77_FLAGS = -fPIC -fdefault-real-8 -fdefault-double-8 -g -O2 -fallow-argument-mismatch
FF90_FLAGS = $(FF77_FLAGS) -std=f2008
FFXX_OPT_FLAGS = -O2
C_FLAGS   = -fPIC -O2 -std=c99

AR       = ar
AR_FLAGS = -rvs

LINKER       = $(FF90)
LINKER_FLAGS =

CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

PYTHON = python3
PYTHON-CONFIG = python3-config
F2PY = f2py
EOF

# 5. 编译
make clean && make

# 6. 安装（需要网络）
pip install -e .

# 如果 HPC 无外网，使用 PYTHONPATH 替代：
# export PYTHONPATH=$(pwd):$PYTHONPATH
```

#### HPC 常见问题

| 问题 | 原因 | 解决方案 |
|------|------|----------|
| `cgns.mod: Unexpected EOF` | `.mod` 文件由 Intel 编译，GCC 不兼容 | 用 GCC 重新编译 CGNS |
| `u_int64_t` 未定义 | `adStack.c` 使用非标准类型 | `#include <stdint.h>`，替换为 `uint64_t` |
| `__intel_sse2_*` 符号未定义 | CGNS 静态库包含 Intel 编译的对象 | 用 GCC 重新编译 CGNS 所有 C 源文件 |
| conda MPI 库冲突 | conda 中的 OpenMPI 5.x 覆盖系统 4.x | 使用编译器绝对路径或 `LD_PRELOAD` |
| `pip install` 失败（无外网） | 无法下载 setuptools | 使用 `PYTHONPATH` 替代 |

### 方式 C：Singularity（支持容器的超算）

```bash
singularity pull adflow-turb.sif docker://mdolab/public:u22-gcc-ompi-stable
```

或使用 Dockerfile 构建自定义镜像：

```dockerfile
FROM mdolab/public:u22-gcc-ompi-stable
COPY patch_adflow_turb.py adflow_turb_ctypes.py /tmp/
RUN python3 /tmp/patch_adflow_turb.py \
    && cd /home/mdolabuser/repos/adflow \
    && make clean && make \
    && pip install -e .
```

## 补丁文件（4 个）

| 文件 | 修改内容 |
|------|----------|
| `src/modules/paramTurb.F90` | SA + SST 常数改为可变；新增 4 个 setter 子程序 |
| `src/initFlow/initializeFlow.F90` | 确保 `referenceState` 中不调用 setter（避免覆盖用户系数） |
| `src/turbulence/SST.F90` | f1 混合函数中硬编码 `0.09`（beta*）替换为 `rSSTBetas` |
| `src/turbulence/turbUtils.F90` | 涡粘计算中硬编码 `0.09`（beta*）替换为 `rSSTBetas` |

## 验证

### HPC A/B 测试（ONERA M6，Paracloud 超算，Job 36920821）

在 ONERA M6 机翼上使用 SA 模型进行三组算例测试（M=0.8395，alpha=3.06 deg）。单个 ADFLOW 求解器实例，通过 ctypes API 在运行间切换系数：

| 算例 | cb1 | CL | CD |
|------|-----|--------|---------|
| A（默认） | 0.1355 | 0.2606478129 | 0.0187979813 |
| B（+48%） | 0.2 | 0.2625755177 | 0.0193669582 |
| C（-41%） | 0.08 | 0.2573790690 | 0.0178077544 |

三组算例产生**不同的** CL/CD。趋势与物理一致：
- cb1 增大 → 湍流产生增强 → 边界层增厚 → CL/CD 增大
- cb1 减小 → 湍流产生减弱 → 边界层减薄 → CL/CD 减小

### Docker NACA 0012 Cp 验证

SA 和 SST 模型分别用 4 组系数（1 组默认 + 3 组随机采样）运行，所有算例均收敛且展现系数依赖的 Cp 分布差异。

详见 `examples/NACA0012/` 目录下的脚本和结果。

## 已知问题

### f2py 模块变量内存复制 Bug（严重）

**不要使用 f2py 模块接口（`pt = solver.adflow.paramturb`）修改系数。** 这是本分支使用 ctypes 的原因。

**症状**：通过 f2py 设置系数后读取值正确，但对 CFD 计算**完全无效**。无论系数如何修改，所有 CL/CD 值精确到机器精度完全相同。

**根因**：f2py 封装 Fortran 共享库时，模块变量数据存在于**两个独立的内存位置**：一个被 f2py Python 封装器访问，另一个被 Fortran 计算子程序（`sa_block`、`blockette` 等）使用。Python 写入副本 A，Fortran 读取副本 B。

**影响编译器**：Intel ifort/ifx 和 GCC gfortran **均受影响**。这不是编译器特定问题。

**验证方法**：在源码中硬编码 `rsaCb1 = 0.5` 后重新编译，CL/CD 产生变化（Job 36920761），证明 Fortran 模块变量本身工作正常。Bug 仅存在于 f2py 的访问方式。

**解决方案**：ctypes API（`adflow_turb_ctypes.py`）使用 `ctypes.c_double.in_dll()` 定位 `libadflow.so` 中的实际 Fortran 符号地址（如 `__paramturb_MOD_rsacb1`），直接写入计算代码读取的同一块内存。已通过 A/B 测试验证（Job 36920788、36920821）。

### CGNS HDF5 格式

如果 CGNS 编译时未链接 HDF5，只能读取 ADF 格式的 CGNS 文件。建议编译 CGNS 时加上 `-DCGNS_ENABLE_HDF5=ON`。

### PETSc 重复初始化

在一个 Python 进程中创建多个 `ADFLOW` 实例会导致 PETSc `MPI_ABORT`。使用单个求解器实例，为不同系数集创建新的 `AeroProblem` 对象。

## 上游版本

基于 [mdolab/adflow](https://github.com/mdolab/adflow) v2.12.1。

## 许可证

与上游一致：GNU Lesser General Public License (LGPL), version 2.1。详见 [LICENSE.md](LICENSE.md)。
