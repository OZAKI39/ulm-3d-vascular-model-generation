#!/usr/bin/env python3
"""一个红细胞在周期 DPD 液体中被带动的入门示例。

基于 Mirheo 官方提交 8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf 的
  tests/doc_scripts/membranes_solvents.py
  tests/fsi/membrane.dp.py
改编。不是原样官方算例，也不是生理参数标定或物理验证结果。
所有量使用示例的模拟单位；不使用真实血管，没有入口/出口或管壁。

运行方式见配套说明。必须使用同一套 Open MPI 启动两个进程：
一个计算进程使用 GPU，另一个负责后处理输出。

Original Mirheo example components:
Copyright (c) 2019 ETH Zurich ( https://www.ethz.ch/en.html )

Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import argparse
import json
import math
import os
from pathlib import Path


def main() -> None:
    # 1. 路径、参数与启动检查。代码文件应放在项目的 py_scripts 下面。
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="一个 RBC + 周期驱动 DPD 液体")
    parser.add_argument("--output", type=Path, required=True,
                        help="本次运行的独立输出目录；不要复用旧结果目录")
    parser.add_argument("--steps", type=int, default=10000,
                        help="计算步数，默认 10000")
    parser.add_argument("--force", type=float, default=0.1,
                        help="液体的驱动力大小，模拟单位；0 表示无外加驱动")
    parser.add_argument("--mesh", type=Path,
                        default=project / "vendor/Mirheo/data/rbc_mesh.off",
                        help="官方 RBC 膜网格路径")
    args = parser.parse_args()

    if args.steps < 100:
        parser.error("本例至少运行 100 步，便于保存核查结果。")
    if not math.isfinite(args.force) or args.force < 0:
        parser.error("--force 必须是有限的非负数。")
    if not args.mesh.is_file():
        parser.error(f"找不到 RBC 网格：{args.mesh}")
    if int(os.environ.get("OMPI_COMM_WORLD_SIZE", "0")) != 2:
        parser.error("请使用 Open MPI 的 mpirun -np 2 启动本脚本，不要直接点击运行按钮。")

    output = args.output.resolve()
    rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
    if rank == 0 and (output / "case_info.json").exists():
        parser.error("该目录已有本例结果；请使用新的 --output，避免覆盖。")
    for directory in (output, output / "ply", output / "h5"):
        directory.mkdir(parents=True, exist_ok=True)

    # 延迟导入：--help 不会启动 Mirheo，也不要求 GPU。
    try:
        import mirheo as mir
    except ImportError as exc:
        raise SystemExit("Mirheo 导入失败，请先完成编译、安装并激活其 .venv。") from exc

    dt = 0.001
    domain = (16.0, 16.0, 16.0)
    number_density = 8.0
    rc = 1.0

    # 2. 创建协调器，再创建一个可变形 RBC 膜。
    u = mir.Mirheo((1, 1, 1), domain,
                   debug_level=2, log_filename=str(output / "log"))
    mesh = mir.ParticleVectors.MembraneMesh(str(args.mesh.resolve()))
    rbc = mir.ParticleVectors.MembraneVector("rbc", mass=1.0, mesh=mesh)
    # 前三个数是放置位置；后四个数是方向四元数，单位四元数表示不旋转。
    ic_rbc = mir.InitialConditions.Membrane([[4.0, 12.0, 8.0, 1.0, 0.0, 0.0, 0.0]])
    u.registerParticleVector(rbc, ic_rbc)

    # 3. 先在盒内生成液体，再按“是否位于膜内”拆分为内、外液体。
    outer = mir.ParticleVectors.ParticleVector("outer", mass=1.0)
    u.registerParticleVector(outer, mir.InitialConditions.Uniform(number_density))
    checker = mir.BelongingCheckers.Mesh("inside_checker")
    u.registerObjectBelongingChecker(checker, rbc)
    inner = u.applyObjectBelongingChecker(
        checker, outer, correct_every=0, inside="inner")
    # inner 已由上面的接口注册；不要再次生成一批内部液体。

    # 4. 连接膜力、DPD 相互作用和防穿膜处理。
    # 以下膜参数来自官方 membranes_solvents.py，不是患者/动物参数。
    membrane_parameters = {
        "x0": 0.457, "ka_tot": 4900.0, "kv_tot": 7500.0,
        "ka": 5000.0, "ks": 0.0444 / 0.000906667,
        "mpow": 2.0, "gammaC": 52.0, "kBT": 0.0,
        "tot_area": 62.2242, "tot_volume": 26.6649,
        "kb": 44.4444, "theta": 6.97,
    }
    membrane = mir.Interactions.MembraneForces(
        "membrane", "wlc", "Kantor", **membrane_parameters)
    u.registerInteraction(membrane)
    u.setInteraction(membrane, rbc, rbc)

    def connect_dpd(name, pv1, pv2, a, gamma):
        interaction = mir.Interactions.Pairwise(
            name, rc, kind="DPD", a=a, gamma=gamma,
            kBT=1.0, power=0.5)
        u.registerInteraction(interaction)
        u.setInteraction(interaction, pv1, pv2)
        return interaction

    # 保留 Python 句柄；这张连接表沿用官方膜—液体教程。
    dpd_handlers = [
        connect_dpd("dpd_outer", outer, outer, 10.0, 10.0),
        connect_dpd("dpd_inner", inner, inner, 10.0, 20.0),
        connect_dpd("dpd_cross", inner, outer, 10.0, 15.0),
        connect_dpd("dpd_outer_membrane", outer, rbc, 0.0, 15.0),
        connect_dpd("dpd_inner_membrane", inner, rbc, 0.0, 15.0),
    ]
    bouncer = mir.Bouncers.Mesh("membrane_bounce", "bounce_maxwell", kBT=0.5)
    u.registerBouncer(bouncer)
    u.setBouncer(bouncer, rbc, outer)
    u.setBouncer(bouncer, rbc, inner)

    # 5. 驱动液体，不直接规定 RBC 的位置或速度。
    # 此版本 direction='x'：y>Ly/2 时为 +Fx，否则为 -Fx。
    fluid_integrator = mir.Integrators.VelocityVerlet_withPeriodicForce(
        "fluid_vv", force=args.force, direction="x")
    membrane_integrator = mir.Integrators.VelocityVerlet("rbc_vv")
    u.registerIntegrator(fluid_integrator)
    u.registerIntegrator(membrane_integrator)
    u.setIntegrator(fluid_integrator, outer)
    u.setIntegrator(fluid_integrator, inner)
    u.setIntegrator(membrane_integrator, rbc)

    # 6. 输出膜形状、中心轨迹、粗粒化平均流场和运行统计。
    u.registerPlugins(mir.Plugins.createStats("stats", every=500))
    u.registerPlugins(mir.Plugins.createParticleChecker("finite_check", check_every=500))
    u.registerPlugins(mir.Plugins.createDumpMesh(
        "mesh_dump", rbc, 100, str(output / "ply") + "/"))
    u.registerPlugins(mir.Plugins.createDumpObjectStats(
        "rbc_stats", rbc, 100, str(output / "rbc_track.csv")))
    u.registerPlugins(mir.Plugins.createDumpAverage(
        "flow_dump", [outer, inner], 10, 1000, (2.0, 2.0, 2.0),
        ["velocities"], str(output / "h5/flow-")))

    if u.isComputeTask():
        info = {
            "example": "single_rbc_periodic_fluid_demo",
            "units": "tutorial simulation units; NOT SI-calibrated",
            "domain": domain, "dt": dt, "steps": args.steps,
            "force": args.force, "number_density": number_density,
            "rbc_mesh": str(args.mesh.resolve()),
            "rbc_placement": [4.0, 12.0, 8.0],
            "force_convention": "+Fx for y>8; -Fx otherwise; solvent only",
            "physical_validation": "NOT_PERFORMED",
        }
        (output / "case_info.json").write_text(
            json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"开始：1 个 RBC，{args.steps} 步；输出目录：{output}", flush=True)

    # 7. 真正推进时间。没有在 Python 中手动移动红细胞。
    u.run(args.steps, dt=dt)
    if u.isComputeTask():
        print("RBC_HELLO_FINISHED：计算循环结束；这不是物理验证通过标志。", flush=True)
        print(f"完整结果请在 mpirun 正常退出后查看：{output}", flush=True)


if __name__ == "__main__":
    main()
