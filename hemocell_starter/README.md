# HemoCell 独立安装与纯流体案例

固定 HemoCell `5a410848bd5c57d5ae1c171112e78eab4a82e650`，配套 Palabos
`05712164d940a42e06afdd705249912fa0c49f14`，不追踪 master。

```bash
bash /home/lzy/projects/hemocell_starter/scripts/setup_hemocell.sh
```

重复安装验证源码版本、归档 SHA256、反向补丁 dry-run 和安装指纹后复用。
不删除源码/Palabos，不改全局环境、默认编译器、CUDA、MPI 或 Mirheo .venv。
首次只编译 `hemocell` 库和独立 `pure_fluid_benchmark` 目标；2 个编译进程。
构建输出：`build/upstream/libhemocell.a` 和 `build/benchmark/pure_fluid_benchmark`。
运行期间安装脚本与两种求解器共享排他锁，避免编译污染计时。

共同调度、物性、原始运行记录与离线页面位于 `../mirheo_starter`；参阅
`../mirheo_starter/test_code/README_solver_benchmark.md`。

许可证：本项目派生 C++ 案例及 HemoCell、Palabos 均为 AGPL-3.0-or-later，
不是 Mirheo 的 MIT 许可。上游头部保留，完整文本在 `cases/pure_fluid_benchmark/COPYING`。
