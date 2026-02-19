"""
NACA0012 SST 系数修改验证: 1 组默认 + 3 组随机采样, 提取上下表面 Cp 并对比。

用法 (Docker 内):
    mpirun --oversubscribe -np 2 python run_naca0012_sst_verify.py

无需编译, 使用 GitHub 预编译的 libadflow.so。
"""
import os
import re
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

comm = MPI.COMM_WORLD
rank = comm.rank

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

# 生成 3 组随机系数 (每组 9 个参数全部从各自范围内独立采样)
np.random.seed(42)
random_sets = []
for i in range(3):
    s = {}
    for name, (default, lo, hi) in SST_PARAMS.items():
        s[name] = np.random.uniform(lo, hi)
    random_sets.append(s)

# 默认系数
default_set = {name: default for name, (default, lo, hi) in SST_PARAMS.items()}

# 全部 4 组: 第 0 组默认, 第 1-3 组随机
all_sets = [("default", default_set)] + [
    (f"random_{i+1}", random_sets[i]) for i in range(3)
]

if rank == 0:
    print("=" * 70)
    print("  NACA0012 SST Coefficient Verification")
    print("  M=0.75, alpha=1.5 deg, 4 cases (1 default + 3 random)")
    print("=" * 70)
    for label, coeffs in all_sets:
        print(f"\n  [{label}]")
        for name, val in coeffs.items():
            default = SST_PARAMS[name][0]
            diff = f"  (delta={val - default:+.5f})" if label != "default" else ""
            print(f"    {name:8s} = {val:.5f}{diff}")

# ---- 公共求解器选项 ----
solverOptions = {
    "gridFile": "/workspace/repo/examples/NACA0012/mesh/n0012.cgns",
    "outputDirectory": "/workspace/repo/examples/NACA0012/output",

    # 物理
    "equationType": "RANS",
    "turbulenceModel": "Menter SST",

    # 求解器 (与官方 tutorial 一致, 仅 turbulenceModel 改为 SST)
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
    "printAllOptions": False,
    "printTiming": False,

    # 输出
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
    "numberSolutions": False,
}

os.makedirs("/workspace/repo/examples/NACA0012/output", exist_ok=True)

# ---- 创建求解器 (只创建一次) ----
solver = ADFLOW(options=solverOptions, comm=comm)
pt = solver.adflow.paramturb

# ---- 确定 z 切片位置 (2D 网格取中间 z 截面) ----
pts_init = solver.getSurfaceCoordinates("wall")
if pts_init is not None and len(pts_init) > 0:
    z_min_local = pts_init[:, 2].min()
    z_max_local = pts_init[:, 2].max()
else:
    z_min_local = 0.0
    z_max_local = 0.0
z_min = comm.allreduce(z_min_local, op=MPI.MIN)
z_max = comm.allreduce(z_max_local, op=MPI.MAX)
z_mid = 0.5 * (z_min + z_max)

if rank == 0:
    print(f"\nMesh z range: [{z_min:.6f}, {z_max:.6f}], slice at z={z_mid:.6f}")

# 添加 z 方向切片 (提取 2D 翼型截面)
solver.addSlices("z", [z_mid])


# ---- 解析 ADflow 切片 .dat 文件 ----
def parse_slice_dat(filepath):
    """解析 ADflow writeSlicesFile 输出的 ASCII Tecplot .dat 文件。"""
    with open(filepath, "r") as f:
        lines = f.readlines()

    # 找变量名行
    var_names = []
    for line in lines:
        if line.strip().upper().startswith("VARIABLES") or line.strip().upper().startswith("\""):
            var_names = re.findall(r'"([^"]+)"', line)
            if var_names:
                break

    if not var_names:
        # 第二行通常是变量定义
        if len(lines) > 1:
            var_names = re.findall(r'"([^"]+)"', lines[1])

    # 解析数据
    data_rows = []
    in_data = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.upper().startswith("ZONE"):
            in_data = True
            continue
        if stripped.upper().startswith("VARIABLES") or stripped.upper().startswith("TITLE"):
            continue
        if stripped.upper().startswith("DATAPACKING") or stripped.upper().startswith("NODES"):
            continue
        if in_data:
            try:
                vals = [float(v) for v in stripped.split()]
                if len(vals) >= 3:
                    data_rows.append(vals)
            except ValueError:
                continue

    if not data_rows:
        return None

    data = np.array(data_rows)
    return {"var_names": var_names, "data": data}


# ---- 逐组求解并提取 Cp ----
results = {}

for case_idx, (label, coeffs) in enumerate(all_sets):
    if rank == 0:
        print(f"\n{'=' * 70}")
        print(f"  Case {case_idx}: {label}")
        print(f"{'=' * 70}")

    # 设置 SST 系数 (9 个参数全部设置)
    pt.setsstconstants(
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

    # 验证系数已写入
    if rank == 0:
        print(f"  Set: sstk={pt.rsstk:.5f}, a1={pt.rssta1:.5f}, "
              f"betas={pt.rsstbetas:.5f}, beta1={pt.rsstbeta1:.5f}, "
              f"sigk2={pt.rsstsigk2:.5f}, beta2={pt.rsstbeta2:.5f}")

    # 设置气动问题
    ap = AeroProblem(
        name=f"n0012_{label}",
        mach=0.75,
        altitude=10000,
        alpha=1.5,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=["cl", "cd"],
    )

    # 求解
    solver(ap)
    funcs = {}
    solver.evalFunctions(ap, funcs)

    cl = funcs.get(f"n0012_{label}_cl", 0.0)
    cd = funcs.get(f"n0012_{label}_cd", 0.0)

    if rank == 0:
        print(f"  Result: Cl = {cl:.6f}, Cd = {cd:.6f}")

    # 切片文件由 solver(ap) 自动写入 (不要显式调用 writeSlicesFile, 会因 famList bug 崩溃)
    # 自动写入路径: /workspace/repo/examples/NACA0012/output/{ap.name}_slices.dat
    if rank == 0:
        slice_file = f"/workspace/repo/examples/NACA0012/output/n0012_{label}_slices.dat"
        if not os.path.exists(slice_file):
            # 尝试找到文件 (可能带不同后缀)
            for fn in os.listdir("/workspace/repo/examples/NACA0012/output"):
                if label in fn and "slice" in fn and fn.endswith(".dat"):
                    slice_file = f"/workspace/repo/examples/NACA0012/output/{fn}"
                    break

        if os.path.exists(slice_file):
            fsize = os.path.getsize(slice_file)
            print(f"  Slice file: {slice_file} ({fsize} bytes)")

            parsed = parse_slice_dat(slice_file)
            if parsed is None:
                print(f"  WARNING: No data in slice file")
                continue

            var_names = parsed["var_names"]
            data = parsed["data"]
            print(f"  Variables: {var_names}")
            print(f"  Data shape: {data.shape}")

            # 查找列索引
            def find_idx(candidates):
                for c in candidates:
                    for k, n in enumerate(var_names):
                        if n.lower().replace(" ", "") == c.lower():
                            return k
                return None

            ix = find_idx(["CoordinateX", "X"])
            iy = find_idx(["CoordinateY", "Y"])
            icp = find_idx(["CoefPressure", "Cp", "cp"])
            ixoc = find_idx(["XoC", "x/c"])

            if icp is None:
                print(f"  WARNING: No Cp column found!")
                continue

            # 用 x/c 列或 CoordinateX 归一化
            if ixoc is not None:
                x = data[:, ixoc]
            elif ix is not None:
                x_raw = data[:, ix]
                x = (x_raw - x_raw.min()) / max(x_raw.max() - x_raw.min(), 1e-10)
            else:
                print(f"  WARNING: No x column found!")
                continue

            cp = data[:, icp]

            # 用 y 坐标区分上下表面
            if iy is not None:
                y = data[:, iy]
                y_mid = np.median(y)
                upper = y >= y_mid
                lower = ~upper
            else:
                # 如果没有 y 列, 用 LE 位置分割
                le_idx = np.argmin(x)
                upper = np.zeros(len(x), dtype=bool)
                lower = np.zeros(len(x), dtype=bool)
                upper[le_idx:] = True
                lower[:le_idx + 1] = True

            x_upper = x[upper]
            cp_upper = cp[upper]
            order_u = np.argsort(x_upper)
            x_upper = x_upper[order_u]
            cp_upper = cp_upper[order_u]

            x_lower = x[lower]
            cp_lower = cp[lower]
            order_l = np.argsort(x_lower)
            x_lower = x_lower[order_l]
            cp_lower = cp_lower[order_l]

            results[label] = {
                "cl": cl, "cd": cd,
                "coeffs": coeffs,
                "x_upper": x_upper, "cp_upper": cp_upper,
                "x_lower": x_lower, "cp_lower": cp_lower,
            }
            print(f"  Surface: {len(x_upper)} upper pts, {len(x_lower)} lower pts")
            print(f"  Cp range: [{cp.min():.4f}, {cp.max():.4f}]")
        else:
            print(f"  WARNING: Slice file not found")

# ---- 保存结果 ----
if rank == 0:
    save_data = {}
    for label, res in results.items():
        prefix = label + "_"
        save_data[prefix + "cl"] = np.float64(res["cl"])
        save_data[prefix + "cd"] = np.float64(res["cd"])
        save_data[prefix + "x_upper"] = res["x_upper"]
        save_data[prefix + "cp_upper"] = res["cp_upper"]
        save_data[prefix + "x_lower"] = res["x_lower"]
        save_data[prefix + "cp_lower"] = res["cp_lower"]
        for pname, pval in res["coeffs"].items():
            save_data[prefix + pname] = np.float64(pval)

    npz_path = "/workspace/repo/examples/NACA0012/output/naca0012_sst_verify.npz"
    np.savez(npz_path, **save_data)
    print(f"\nSaved: {npz_path}")

    # ---- 绘图 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = ["k", "C0", "C1", "C2"]
    lstyles = ["-", "--", "-.", ":"]

    for i, (label, res) in enumerate(results.items()):
        c = colors[i]
        ls = lstyles[i]
        lw = 2.0 if i == 0 else 1.2
        tag = f"{label} (Cl={res['cl']:.4f}, Cd={res['cd']:.5f})"

        ax.plot(res["x_upper"], res["cp_upper"], color=c, linestyle=ls,
                linewidth=lw, label=f"{tag} upper")
        ax.plot(res["x_lower"], res["cp_lower"], color=c, linestyle=ls,
                linewidth=lw, alpha=0.6)

    ax.set_xlabel("$x/c$", fontsize=12)
    ax.set_ylabel("$C_p$", fontsize=12)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    ax.set_title(r"NACA 0012, $M=0.75$, $\alpha=1.5^\circ$ — SST Coefficient Sensitivity",
                 fontsize=13)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()

    fig_path = "/workspace/repo/examples/NACA0012/output/naca0012_sst_cp_compare.png"
    fig.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")

    # 打印对比表
    print(f"\n{'=' * 70}")
    print(f"  Summary")
    print(f"{'=' * 70}")
    print(f"  {'Case':<12s} {'Cl':>10s} {'Cd':>10s}")
    print(f"  {'-'*12} {'-'*10} {'-'*10}")
    for label, res in results.items():
        print(f"  {label:<12s} {res['cl']:10.6f} {res['cd']:10.6f}")

    print("\nDONE")
