"""
NACA0012 SST 系数修改验证 (HPC ctypes 版)
==========================================
1 组默认 + 3 组随机采样 (9 个参数), 提取 CL/CD 并对比。
使用 ctypes API 绕过 f2py 模块变量内存复制 bug。

依赖: 将 adflow_turb_ctypes.py 放在本脚本同目录下，或加入 PYTHONPATH。

用法 (HPC):
    mpirun -np 4 python validate_naca0012_sst.py
"""
import os
import sys
import re
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

# adflow_turb_ctypes.py 应与本脚本在同一目录，或已在 PYTHONPATH 中
from adflow_turb_ctypes import (
    set_sst_constants, set_sst_defaults, get_sst_constants,
)

comm = MPI.COMM_WORLD
rank = comm.rank

# ---- 输出目录 ----
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output_sst")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- 网格路径 (HPC) ----
MESH_FILE = os.environ.get(
    "NACA0012_MESH",
    "/public3/home/m6s001303/HPC_BM/NACA0012/mesh/n0012.cgns",
)

# ---- SST 校准参数范围 (9 个参数全部独立采样) ----
SST_PARAMS = {
    #  name        default   lo      hi
    "sstk":      (0.41,    0.35,   0.45),
    "a1":        (0.31,    0.25,   0.40),
    "betas":     (0.09,    0.075,  0.11),
    "sigk1":     (0.85,    0.60,   1.20),
    "sigw1":     (0.50,    0.30,   0.70),
    "beta1":     (0.075,   0.05,   0.10),
    "sigk2":     (1.00,    0.70,   1.30),
    "sigw2":     (0.856,   0.60,   1.10),
    "beta2":     (0.0828,  0.06,   0.11),
}

# 生成 3 组随机系数
np.random.seed(42)
random_sets = []
for i in range(3):
    s = {}
    for name, (default, lo, hi) in SST_PARAMS.items():
        s[name] = np.random.uniform(lo, hi)
    random_sets.append(s)

# 默认系数
default_set = {name: default for name, (default, lo, hi) in SST_PARAMS.items()}

# 全部 4 组
all_sets = [("default", default_set)] + [
    (f"random_{i+1}", random_sets[i]) for i in range(3)
]

if rank == 0:
    print("=" * 70)
    print("  NACA0012 SST Coefficient Verification (ctypes API, HPC)")
    print("  M=0.75, alpha=1.5 deg, 4 cases (1 default + 3 random)")
    print("=" * 70)
    print(f"  Mesh: {MESH_FILE}")
    print(f"  Output: {OUTPUT_DIR}")
    for label, coeffs in all_sets:
        print(f"\n  [{label}]")
        for name, val in coeffs.items():
            default = SST_PARAMS[name][0]
            diff = f"  (delta={val - default:+.5f})" if label != "default" else ""
            print(f"    {name:8s} = {val:.5f}{diff}")

# ---- 求解器选项 ----
solverOptions = {
    "gridFile": MESH_FILE,
    "outputDirectory": OUTPUT_DIR,
    "equationType": "RANS",
    "turbulenceModel": "Menter SST",
    "MGCycle": "sg",
    "nSubiterTurb": 10,
    "useANKSolver": True,
    "useNKSolver": True,
    "NKSwitchTol": 1e-4,
    "L2Convergence": 1e-15,
    "nCycles": 1000,
    "monitorVariables": ["resrho", "cl", "cd"],
    "printIterations": True,
    "printAllOptions": False,
    "printTiming": False,
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
    "numberSolutions": False,
}

# ---- 创建求解器 ----
solver = ADFLOW(options=solverOptions, comm=comm)

# ---- z 切片 ----
pts = solver.getSurfaceCoordinates("wall")
if pts is not None and len(pts) > 0:
    z_min_l, z_max_l = pts[:, 2].min(), pts[:, 2].max()
else:
    z_min_l, z_max_l = 0.0, 0.0
z_min = comm.allreduce(z_min_l, op=MPI.MIN)
z_max = comm.allreduce(z_max_l, op=MPI.MAX)
z_mid = 0.5 * (z_min + z_max)
if rank == 0:
    print(f"\nMesh z: [{z_min:.6f}, {z_max:.6f}], slice at z={z_mid:.6f}")
solver.addSlices("z", [z_mid])


# ---- 解析切片 .dat ----
def parse_slice_dat(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    var_names = []
    for line in lines:
        if "VARIABLES" in line.upper():
            var_names = re.findall(r'"([^"]+)"', line)
            if var_names:
                break
    data_rows = []
    in_data = False
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.upper().startswith("ZONE"):
            in_data = True
            continue
        if s.upper().startswith(("VARIABLES", "TITLE", "DATAPACKING")):
            continue
        if in_data:
            try:
                vals = [float(v) for v in s.split()]
                if len(vals) >= 3:
                    data_rows.append(vals)
            except ValueError:
                continue
    if not data_rows:
        return None
    return {"var_names": var_names, "data": np.array(data_rows)}


# ---- 逐组求解 ----
results = {}

for case_idx, (label, coeffs) in enumerate(all_sets):
    if rank == 0:
        print(f"\n{'=' * 70}")
        print(f"  Case {case_idx}: {label}")
        print(f"{'=' * 70}")

    # 通过 ctypes 设置 SST 系数
    set_sst_constants(
        sstk=coeffs["sstk"],
        a1=coeffs["a1"],
        betas=coeffs["betas"],
        sigk1=coeffs["sigk1"],
        sigw1=coeffs["sigw1"],
        beta1=coeffs["beta1"],
        sigk2=coeffs["sigk2"],
        sigw2=coeffs["sigw2"],
        beta2=coeffs["beta2"],
    )

    # 验证
    if rank == 0:
        c = get_sst_constants()
        print(f"  ctypes set: sstk={c['sstk']:.5f}, a1={c['a1']:.5f}, "
              f"betas={c['betas']:.5f}")
        print(f"  ctypes set: sigk1={c['sigk1']:.5f}, sigw1={c['sigw1']:.5f}, "
              f"beta1={c['beta1']:.5f}")
        print(f"  ctypes set: sigk2={c['sigk2']:.5f}, sigw2={c['sigw2']:.5f}, "
              f"beta2={c['beta2']:.5f}")

    ap = AeroProblem(
        name=f"n0012_sst_{label}",
        mach=0.75, altitude=10000, alpha=1.5,
        areaRef=1.0, chordRef=1.0,
        evalFuncs=["cl", "cd"],
    )

    solver(ap)
    funcs = {}
    solver.evalFunctions(ap, funcs)

    cl = funcs.get(f"n0012_sst_{label}_cl", 0.0)
    cd = funcs.get(f"n0012_sst_{label}_cd", 0.0)

    if rank == 0:
        print(f"  Result: CL = {cl:.10f}, CD = {cd:.10f}")

    # 切片文件
    if rank == 0:
        slice_file = os.path.join(OUTPUT_DIR, f"n0012_sst_{label}_slices.dat")
        if not os.path.exists(slice_file):
            for fn in os.listdir(OUTPUT_DIR):
                if label in fn and "slice" in fn and fn.endswith(".dat"):
                    slice_file = os.path.join(OUTPUT_DIR, fn)
                    break

        if os.path.exists(slice_file):
            fsize = os.path.getsize(slice_file)
            print(f"  Slice file: {slice_file} ({fsize} bytes)")

            parsed = parse_slice_dat(slice_file)
            if parsed is None:
                print(f"  WARNING: No data in slice file")
                results[label] = {"cl": cl, "cd": cd, "coeffs": coeffs}
                continue

            var_names = parsed["var_names"]
            data = parsed["data"]

            def find_idx(candidates):
                for c_name in candidates:
                    for k, n in enumerate(var_names):
                        if n.lower().replace(" ", "") == c_name.lower():
                            return k
                return None

            ix = find_idx(["XoC", "CoordinateX", "X"])
            iy = find_idx(["YoC", "CoordinateY", "Y"])
            icp = find_idx(["CoefPressure", "Cp"])

            if icp is not None and ix is not None:
                x = data[:, ix]
                cp = data[:, icp]
                y = data[:, iy] if iy is not None else None

                if y is not None:
                    upper = y >= 0
                    lower = y < 0
                else:
                    le_idx = np.argmin(x)
                    upper = np.arange(len(x)) <= le_idx
                    lower = np.arange(len(x)) >= le_idx

                x_up, cp_up = x[upper], cp[upper]
                order_u = np.argsort(x_up)
                x_up, cp_up = x_up[order_u], cp_up[order_u]

                x_lo, cp_lo = x[lower], cp[lower]
                order_l = np.argsort(x_lo)
                x_lo, cp_lo = x_lo[order_l], cp_lo[order_l]

                results[label] = {
                    "cl": cl, "cd": cd, "coeffs": coeffs,
                    "x_upper": x_up, "cp_upper": cp_up,
                    "x_lower": x_lo, "cp_lower": cp_lo,
                }
                print(f"  Surface: {len(x_up)} upper, {len(x_lo)} lower")
                print(f"  Cp range: [{cp.min():.4f}, {cp.max():.4f}]")
            else:
                print(f"  WARNING: Missing Cp or x column")
                results[label] = {"cl": cl, "cd": cd, "coeffs": coeffs}
        else:
            print(f"  WARNING: Slice file not found")
            results[label] = {"cl": cl, "cd": cd, "coeffs": coeffs}

# ---- 保存 & 汇总 ----
if rank == 0:
    save_data = {}
    for label, res in results.items():
        p = label + "_"
        save_data[p + "cl"] = np.float64(res["cl"])
        save_data[p + "cd"] = np.float64(res["cd"])
        if "x_upper" in res:
            save_data[p + "x_upper"] = res["x_upper"]
            save_data[p + "cp_upper"] = res["cp_upper"]
            save_data[p + "x_lower"] = res["x_lower"]
            save_data[p + "cp_lower"] = res["cp_lower"]
        for pname, pval in res["coeffs"].items():
            save_data[p + pname] = np.float64(pval)

    npz_path = os.path.join(OUTPUT_DIR, "naca0012_sst_verify.npz")
    np.savez(npz_path, **save_data)
    print(f"\nSaved: {npz_path}")

    # 对比表
    print(f"\n{'=' * 70}")
    print(f"  Summary")
    print(f"{'=' * 70}")
    print(f"  {'Case':<12s} {'CL':>12s} {'CD':>12s}")
    print(f"  {'-'*12} {'-'*12} {'-'*12}")
    for label, res in results.items():
        print(f"  {label:<12s} {res['cl']:12.8f} {res['cd']:12.8f}")

    # 验证
    cls = [res["cl"] for res in results.values()]
    if len(set(f"{v:.10f}" for v in cls)) == 1:
        print("\n  !!! WARNING: ALL CL VALUES IDENTICAL - ctypes may not be working !!!")
    else:
        print(f"\n  >>> ALL {len(cls)} CASES DIFFERENT - ctypes API VERIFIED <<<")

    print("\nDONE")
