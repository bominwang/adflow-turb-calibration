# 离线依赖包

本目录包含 ADflow 编译所需的全部依赖源码，用于无法联网的 HPC 环境。

## 包含的源码包

| 文件 | 版本 | 说明 |
|------|------|------|
| `petsc-3.15.5.tar.gz` | 3.15.5 | 并行线性代数求解器 |
| `CGNS-4.2.0.tar.gz` | 4.2.0 | CFD 网格文件格式库 |
| `superlu_dist-6.4.0.tar.gz` | 6.4.0 | 分布式稀疏直接求解器（PETSc NK 求解器加速） |

## 使用方法

```bash
# 编译全部依赖（PETSc + CGNS）
bash hpc/build_deps.sh

# 或指定安装路径
bash hpc/build_deps.sh /path/to/install

# 然后编译 ADflow
export PETSC_DIR=./deps/install/petsc
export PETSC_ARCH=""
export CGNS_HOME=./deps/install/cgns
bash hpc/deploy.sh
```

## 注意事项

- SuperLU_DIST 会被 PETSc 自动编译和集成
- CGNS 默认编译为静态库（不依赖 HDF5）；如需 HDF5 格式支持，需自行安装 HDF5 并修改 cmake 参数
- 如果 HPC 上已有 PETSc 和 CGNS，可以直接设置环境变量，跳过此步骤
