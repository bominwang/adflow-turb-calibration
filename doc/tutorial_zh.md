# ADflow 湍流模型系数修改教程

本教程面向初学者，介绍如何使用修改版 ADflow 进行 SA/SST 湍流模型闭合系数的运行时修改。适用于贝叶斯校准、不确定性量化等研究场景。

## 目录

- [背景知识](#背景知识)
- [环境准备](#环境准备)
- [快速开始](#快速开始)
- [SA 模型系数修改](#sa-模型系数修改)
- [SST 模型系数修改](#sst-模型系数修改)
- [完整算例：NACA 0012](#完整算例naca-0012)
- [在 HPC 上运行](#在-hpc-上运行)
- [常见问题](#常见问题)
- [API 参考](#api-参考)

---

## 背景知识

### 为什么要修改湍流模型系数

RANS（雷诺平均 Navier-Stokes）方程中的湍流模型（如 SA、SST）包含若干**闭合系数**。这些系数是通过基准实验确定的经验值，但在不同流动条件下可能不是最优的。

通过贝叶斯校准，可以：
- 量化闭合系数的不确定性
- 根据实验数据更新系数的后验分布
- 评估模型预测的置信区间

### SA 和 SST 模型简介

**Spalart-Allmaras (SA)** 是单方程湍流模型，求解湍流粘性的输运方程。包含 9 个可校准参数。

**Menter SST (k-omega SST)** 是两方程湍流模型，结合了 k-epsilon 在远场和 k-omega 在近壁面的优势。包含 9 个独立参数。

### 原理概述

原版 ADflow 将这些系数声明为 Fortran `parameter`（编译期常量）。我们的修改：

1. 将系数从 `parameter` 改为普通模块变量（运行时可写）
2. 通过 Python ctypes 直接写入 Fortran 共享库中的符号地址
3. 修改后的系数**立即生效**于下一次 CFD 求解

```
Python (ctypes)  →  直接写入 libadflow.so 中的内存地址  →  Fortran 计算代码读取
```

> **为什么不用 f2py？** f2py 封装 Fortran 模块时会产生**内存复制**：Python 端和 Fortran 端各自持有一份模块变量副本。Python 写入的是 f2py 的副本，Fortran 计算代码读取的是自己的副本。结果是修改看似成功（读回值正确），但对 CFD 计算**完全无效**。这个 bug 同时影响 Intel ifort 和 GCC gfortran。

---

## 环境准备

### 方式一：Docker（推荐初学者）

```bash
# 1. 拉取镜像
docker pull mdolab/public:u22-gcc-ompi-stable

# 2. 启动容器
docker run -it --name adflow-turb mdolab/public:u22-gcc-ompi-stable bash

# 3. 在容器内克隆仓库并补丁
git clone https://github.com/bominwang/adflow-turb-calibration.git /tmp/repo
python3 /tmp/repo/patch_adflow_turb.py

# 4. 编译（约 5-10 分钟）
cd /home/mdolabuser/repos/adflow
make clean && make
pip install -e .

# 5. 复制 ctypes API 到工作目录
cp /tmp/repo/adflow_turb_ctypes.py /home/mdolabuser/
```

编译成功的标志：
```
Testing if module libadflow can be imported...
Module libadflow was successfully imported
```

### 方式二：HPC 原生编译

参见 [README.md](README.md) 中 "方式 B：HPC 原生编译" 章节，或使用 `hpc/deploy.sh` 一键部署脚本。

---

## 快速开始

### 第一个例子：修改 SA 的 cb1 系数

```python
from adflow import ADFLOW
from adflow_turb_ctypes import set_sa_constants, get_sa_constants, set_sa_defaults

# 1. 创建求解器（只创建一次）
options = {
    "gridFile": "your_mesh.cgns",
    "outputDirectory": "./output",
    "equationType": "RANS",
    # SA 是 ADflow 默认湍流模型，不需要额外指定
}
solver = ADFLOW(options=options)

# 2. 查看当前系数（默认值）
print(get_sa_constants())
# 输出: {'cb1': 0.1355, 'cb2': 0.622, 'sigma': 0.6667, ...}

# 3. 修改 cb1（其余保持默认）
set_sa_constants(
    cb1=0.2,           # 修改这个
    cb2=0.622,         # 以下保持默认
    sigma=2.0/3.0,
    kappa=0.41,
    cv1=7.1,
    cw2=0.3,
    cw3=2.0,
    ct3=1.2,
    ct4=0.5,
)

# 4. 验证修改生效
c = get_sa_constants()
print(f"cb1 = {c['cb1']}")  # 应该输出 0.2
print(f"cw1 = {c['cw1']}")  # cw1 会自动重算

# 5. 求解
from baseclasses import AeroProblem
ap = AeroProblem(name="test", mach=0.8, alpha=2.0, altitude=10000,
                 areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd"])
solver(ap)

funcs = {}
solver.evalFunctions(ap, funcs)
print(f"CL = {funcs['test_cl']:.6f}")
print(f"CD = {funcs['test_cd']:.6f}")

# 6. 恢复默认
set_sa_defaults()
```

### 核心原则

1. **求解器只创建一次**。创建多个 `ADFLOW` 实例会导致 PETSc 崩溃。
2. **每次求解前设置系数**。ADflow 切换 AeroProblem 时可能重置内部状态。
3. **adflow_turb_ctypes.py 和脚本放在同一目录**。Python 自动将脚本所在目录加入搜索路径。

---

## SA 模型系数修改

### 9 个可校准参数

```python
set_sa_constants(
    cb1=0.1355,       # 产生项系数
    cb2=0.622,        # 扩散项系数
    sigma=2.0/3.0,    # 扩散比（注意：Fortran 中存为 rsaCb3）
    kappa=0.41,       # von Karman 常数
    cv1=7.1,          # 壁面阻尼系数
    cw2=0.3,          # 耗散项系数
    cw3=2.0,          # 耗散项系数
    ct3=1.2,          # 转捩函数系数
    ct4=0.5,          # 转捩函数系数
)
```

### 派生量 cw1

`cw1` 是根据其他系数自动计算的：

```
cw1 = cb1 / kappa^2 + (1 + cb2) / sigma
```

默认值：`cw1 = 0.1355 / 0.41^2 + 1.622 / 0.6667 = 3.2391`

每次调用 `set_sa_constants()` 后，cw1 会自动重算。

### 常见校准策略

通常只校准 7 个参数（不含 ct3/ct4），因为 ct3 和 ct4 主要影响转捩行为，对全湍流计算影响极小：

```python
# 7 参数校准（ct3/ct4 保持默认）
set_sa_constants(
    cb1=new_cb1, cb2=new_cb2, sigma=new_sigma,
    kappa=new_kappa, cv1=new_cv1, cw2=new_cw2, cw3=new_cw3,
    ct3=1.2,    # 默认
    ct4=0.5,    # 默认
)
```

### 参数物理范围（建议）

| 参数 | 默认值 | 建议范围 | 物理约束 |
|------|--------|----------|----------|
| cb1 | 0.1355 | [0.12, 0.15] | 正值 |
| cb2 | 0.622 | [0.55, 0.75] | 正值 |
| sigma | 2/3 | [0.5, 1.5] | 正值 |
| kappa | 0.41 | [0.36, 0.44] | 正值，对数律中的 von Karman 常数 |
| cv1 | 7.1 | [6.5, 8.0] | 正值 |
| cw2 | 0.3 | [0.05, 0.4] | [0, 1) |
| cw3 | 2.0 | [1.5, 3.0] | 正值 |

---

## SST 模型系数修改

### 9 个独立参数

```python
from adflow_turb_ctypes import set_sst_constants, get_sst_constants, set_sst_defaults

set_sst_constants(
    sstk=0.41,        # kappa (von Karman 常数)
    a1=0.31,          # 限制器常数
    betas=0.09,       # beta* (耗散系数)
    sigk1=0.85,       # k 方程 set 1 扩散系数
    sigw1=0.5,        # omega 方程 set 1 扩散系数
    beta1=0.075,      # set 1 耗散系数
    sigk2=1.0,        # k 方程 set 2 扩散系数
    sigw2=0.856,      # omega 方程 set 2 扩散系数
    beta2=0.0828,     # set 2 耗散系数
)
```

### 使用 SST 模型

ADflow 默认使用 SA。要切换到 SST，需要在求解器选项中指定：

```python
options = {
    "gridFile": "your_mesh.cgns",
    "equationType": "RANS",
    "turbulenceModel": "Menter SST",   # 关键：指定 SST
    # ... 其他选项
}
solver = ADFLOW(options=options)
```

---

## 完整算例：NACA 0012

### 问题描述

NACA 0012 翼型，M=0.75，迎角 1.5 度。比较默认系数和修改后系数的 CL/CD 差异。

### 求解器配置

```python
import os
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from adflow_turb_ctypes import set_sa_constants, get_sa_constants, set_sa_defaults

comm = MPI.COMM_WORLD

options = {
    "gridFile": "n0012.cgns",
    "outputDirectory": "./output",
    "equationType": "RANS",

    # 求解器设置
    "MGCycle": "sg",
    "nSubiterTurb": 10,
    "useANKSolver": True,
    "useNKSolver": True,
    "NKSwitchTol": 1e-4,
    "L2Convergence": 1e-15,
    "nCycles": 1000,

    # 监控
    "monitorVariables": ["resrho", "cl", "cd"],
    "printIterations": True,

    # 输出控制
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
}

solver = ADFLOW(options=options, comm=comm)
```

### 多组系数对比

```python
# 定义多组系数
cases = {
    "default": {"cb1": 0.1355},
    "cb1_high": {"cb1": 0.2},
    "cb1_low": {"cb1": 0.08},
}

results = {}
for name, params in cases.items():
    # 先恢复默认，再修改指定参数
    set_sa_defaults()
    defaults = get_sa_constants()

    set_sa_constants(
        cb1=params.get("cb1", defaults["cb1"]),
        cb2=params.get("cb2", defaults["cb2"]),
        sigma=params.get("sigma", defaults["sigma"]),
        kappa=params.get("kappa", defaults["kappa"]),
        cv1=params.get("cv1", defaults["cv1"]),
        cw2=params.get("cw2", defaults["cw2"]),
        cw3=params.get("cw3", defaults["cw3"]),
        ct3=params.get("ct3", defaults["ct3"]),
        ct4=params.get("ct4", defaults["ct4"]),
    )

    ap = AeroProblem(
        name=f"naca0012_{name}",
        mach=0.75, altitude=10000, alpha=1.5,
        areaRef=1.0, chordRef=1.0,
        evalFuncs=["cl", "cd"],
    )

    solver(ap)
    funcs = {}
    solver.evalFunctions(ap, funcs)

    results[name] = {
        "cl": funcs[f"naca0012_{name}_cl"],
        "cd": funcs[f"naca0012_{name}_cd"],
    }

# 打印结果
if comm.rank == 0:
    print(f"{'Case':<12} {'CL':>10} {'CD':>10}")
    print("-" * 34)
    for name, r in results.items():
        print(f"{name:<12} {r['cl']:10.6f} {r['cd']:10.6f}")
```

### 提取 Cp 分布

```python
# 添加 z 方向切片（2D 翼型取中间截面）
pts = solver.getSurfaceCoordinates("wall")
z_mid = 0.5 * (pts[:, 2].min() + pts[:, 2].max())
solver.addSlices("z", [z_mid])

# 求解后自动写入切片文件
solver(ap)
# 切片文件: output/naca0012_default_slices.dat
```

> **注意**：不要手动调用 `solver.writeSlicesFile()`，它有一个 famList 参数 bug。切片会在 `solver(ap)` 时自动写入。

---

## 在 HPC 上运行

### Slurm 作业脚本模板

```bash
#!/bin/bash
#SBATCH -J adflow_sa        # 作业名
#SBATCH -p amd_512           # 分区（根据 HPC 修改）
#SBATCH -N 1                 # 节点数
#SBATCH -n 4                 # MPI 进程数
#SBATCH -o result.out        # 标准输出
#SBATCH -e result.err        # 标准错误
#SBATCH -t 00:30:00          # 最大运行时间

# 加载环境（根据 HPC 修改）
source /path/to/env_gcc.sh
export PYTHONPATH=/path/to/adflow-gcc:$PYTHONPATH

# 运行
cd /path/to/your/scripts
mpirun -np 4 python your_script.py
```

### 文件部署

将以下文件放在同一目录：
```
your_working_dir/
├── adflow_turb_ctypes.py    # ctypes API（从仓库复制）
├── your_script.py           # 你的计算脚本
├── your_mesh.cgns           # 网格文件
└── job.slurm                # 作业提交脚本
```

### 提交和监控

```bash
# 提交作业
sbatch job.slurm

# 查看队列
squeue -u $USER

# 查看输出（作业运行中）
tail -f result.out

# 取消作业
scancel <job_id>
```

---

## 常见问题

### Q: 修改系数后 CL/CD 没有变化？

检查是否使用了 ctypes API (`adflow_turb_ctypes.py`)。**不要**使用 f2py 接口 (`solver.adflow.paramturb`)，它有内存复制 bug，修改不会生效。

### Q: 报错 `ModuleNotFoundError: No module named 'adflow_turb_ctypes'`？

确保 `adflow_turb_ctypes.py` 与你的脚本在**同一目录**下，或者已加入 `PYTHONPATH`。

### Q: 报错 PETSc MPI_ABORT？

不要创建多个 `ADFLOW` 实例。使用单个求解器，创建不同的 `AeroProblem` 来切换工况。

### Q: 可以只修改某一个系数吗？

`set_sa_constants()` 和 `set_sst_constants()` 要求传入全部参数。推荐做法：

```python
# 先获取当前值
c = get_sa_constants()
# 只修改你想改的
set_sa_constants(
    cb1=0.2,              # 修改这个
    cb2=c['cb2'],         # 其余保持不变
    sigma=c['sigma'],
    kappa=c['kappa'],
    cv1=c['cv1'],
    cw2=c['cw2'],
    cw3=c['cw3'],
    ct3=c['ct3'],
    ct4=c['ct4'],
)
```

### Q: 编译时报 `Type mismatch` 警告？

这是 ADflow 上游 CGNS 模块的警告，不影响计算结果。只要最终显示 `Module libadflow was successfully imported` 就表示编译成功。

### Q: `cgns.mod: Unexpected EOF` 错误？

CGNS 的 `.mod` 文件是由不同编译器生成的（比如 Intel 编译的 `.mod` 不能被 GCC 读取）。需要用相同的编译器重新编译 CGNS。

### Q: 怎么判断 ctypes API 是否正常工作？

运行简单的读写测试：
```python
from adflow_turb_ctypes import set_sa_constants, get_sa_constants, set_sa_defaults

set_sa_defaults()
c1 = get_sa_constants()
print(f"默认 cb1 = {c1['cb1']}")   # 应为 0.1355

set_sa_constants(0.5, 0.622, 2/3, 0.41, 7.1, 0.3, 2.0, 1.2, 0.5)
c2 = get_sa_constants()
print(f"修改后 cb1 = {c2['cb1']}")  # 应为 0.5

set_sa_defaults()
c3 = get_sa_constants()
print(f"重置后 cb1 = {c3['cb1']}")  # 应为 0.1355
```

---

## API 参考

### SA 模型

```python
from adflow_turb_ctypes import (
    set_sa_constants,    # 设置 9 个 SA 参数（cw1 自动重算）
    set_sa_defaults,     # 恢复 SA 默认值
    get_sa_constants,    # 读取当前 SA 系数（返回 dict）
)
```

**`set_sa_constants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)`**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| cb1 | float | 0.1355 | 产生项系数 |
| cb2 | float | 0.622 | 扩散项系数 |
| sigma | float | 2/3 | 扩散比 |
| kappa | float | 0.41 | von Karman 常数 |
| cv1 | float | 7.1 | 壁面阻尼系数 |
| cw2 | float | 0.3 | 耗散项系数 |
| cw3 | float | 2.0 | 耗散项系数 |
| ct3 | float | 1.2 | 转捩函数系数 |
| ct4 | float | 0.5 | 转捩函数系数 |

**`get_sa_constants()`** → `dict`

返回包含 10 个键的字典：9 个输入参数 + 派生量 `cw1`。

**`set_sa_defaults()`**

将所有 SA 系数恢复为文献标准默认值。

### SST 模型

```python
from adflow_turb_ctypes import (
    set_sst_constants,   # 设置 9 个 SST 参数
    set_sst_defaults,    # 恢复 SST 默认值
    get_sst_constants,   # 读取当前 SST 系数（返回 dict）
)
```

**`set_sst_constants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)`**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| sstk | float | 0.41 | von Karman 常数 |
| a1 | float | 0.31 | 限制器常数 |
| betas | float | 0.09 | beta* 耗散系数 |
| sigk1 | float | 0.85 | set 1 k 扩散系数 |
| sigw1 | float | 0.5 | set 1 omega 扩散系数 |
| beta1 | float | 0.075 | set 1 耗散系数 |
| sigk2 | float | 1.0 | set 2 k 扩散系数 |
| sigw2 | float | 0.856 | set 2 omega 扩散系数 |
| beta2 | float | 0.0828 | set 2 耗散系数 |

---

## 仓库结构

```
adflow-turb-calibration/
├── adflow_turb_ctypes.py        # ★ ctypes Python API（核心文件）
├── patch_adflow_turb.py         # 补丁脚本（修改 4 个 Fortran 文件）
├── test_turb_coefficients.py    # 系数读写验证脚本
├── README.md                    # 项目说明
├── src/                         # ADflow 源码
├── examples/
│   └── NACA0012/
│       ├── mesh/                # 网格文件
│       ├── run_naca0012_sa_verify.py
│       └── validation/          # HPC 验证脚本
├── hpc/
│   ├── deploy.sh               # 一键部署脚本
│   ├── build_deps.sh           # 离线依赖编译脚本
│   └── env_gcc_template.sh     # 环境模板
└── deps/                        # 离线依赖源码（offline 分支）
```

---

*最后更新: 2026-02-21*
