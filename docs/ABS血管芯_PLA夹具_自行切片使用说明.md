# ABS 血管芯、PLA 夹具：自行切片

**双上端口芯最新姿态选择：** 优先查看 `deliverables/Straight_2Top_P1S04_ABS_orientation_options_unsliced.zip`，内有少支撑的约 45° 斜放方案，以及主直管平行打印板、离板约 1.20 mm 的方案，两者择一打印。仍不切片，详见 [姿态选择说明](Straight_2Top少支撑与水平离板姿态.md)。下方竖直版和原六份包保留供追溯，其余 PLA 文件不受影响。

**2026-10-15 前一版精细打印方案：** `deliverables/Straight_2Top_P1S04_ABS_fine_vertical_unsliced.zip` 采用主直管竖直、0.08 mm 层高、低速及无筏层设置，保留用于比较；详见 [精细打印说明](Straight_2Top血管芯精细打印.md)。下述六份工程包保留原版姿态，其他五份仍照常使用。

2026-10-12 当前交付方式。三版模型均采用 **Bambu P1S / 0.4 mm 喷嘴**；只生成模型和可编辑、未切片的 3MF，不调用切片引擎、不生成打印刀路、不提供预计打印时间。

2026-10-12 名称兼容性修复：打印盘名称改用连字符，移除了 Bambu Studio 禁止的 `|`。导出时检查项目、打印盘和零件名称，包括 XML 解码后的非法字符；此前已下载的工程请用本次更新包替换。模型尺寸、材料与排版不变。

2026-10-12 图形界面兼容性修复：补齐工程级 `filament_colour`、`filament_ids`、`filament_self_index` 等字段。预设中的 Standard / High Flow 是同一耗材的挤出变体，需要映射到同一个料盘；缺失此映射时，命令行仍可读出网格，但图形界面可能先报“无效配置”，再报“不包含任何几何数据”。本次以实际 GUI 打开验证为准，未执行切片。最新六份工程汇总包为 `deliverables/P1S_04_ABS_PLA_unsliced_3MF_GUI_verified.zip`；此前的 `name_fixed.zip` 已归档。

| 组件 | 材料 |
| --- | --- |
| 血管芯，包括竖直端口与定位肩台 | ABS |
| 无底四壁盒、外伸安装脚、燕尾座、两侧磁铁槽 | PLA |
| 独立燕尾支座 | PLA |
| 横向锁键 | PLA |

每套模型提供两个分开的单材料工程：`*_core_ABS_project.3mf` 只有血管芯，`*_fixture_PLA_project.3mf` 包含盒体、支座和锁键。分别打开工程、确认对应料盘与打印参数后自行切片，不需要在同一次打印中混用 ABS 和 PLA。原装配模型仍完整保留，生产件数不变：Y 形 8 件，两种直管各 6 件。

## 模型与尺寸

本轮改变材料分配和交付方式，血管芯、框架、磁铁槽、支座、锁键的实体形状不改变。血管芯直径 1.60 mm，距外部浇注底面的最薄 PDMS 层为 1.50 mm，盒壁 2.00 mm，盒内无打印底板。磁铁槽内径 6.40 mm、深 2.00 mm、底厚 1.00 mm。盒体顶高 10.40 mm，支座顶高 17.80 mm，血管芯顶高 18.40 mm；托臂与盒沿仍相隔 3.00 mm。

腔体原有的 CAD 宽/高下限 9.60 mm 保留，不再将 ABS 的 2% 收缩预算作为 PLA 的收缩参数。配置使用 `cad_dimensional_allowance` 明确保存 9.60 mm 既有 CAD 下限和 0.30 mm 未实物校准的尺寸误差预算，没有假设 PLA 的实际收缩率。工程默认保持 100% 缩放和零额外孔/轮廓补偿；实际打印尺寸仍需自行校准。

3dm 为独立命名、分层的封闭网格，毫米单位，并在组件属性中记录 ABS 或 PLA。血管芯另有解析 STEP；各组件都有 STL。装配 3dm 用于查看部件关系，分件排版 3dm 用于查看打印姿态；实际切片请优先打开按材料拆分的 3MF。

## 当前输出

以下路径相对于 `/home/lzy/projects/vascular_printing`；服务器项目根为 `/workspace/vascular_printing`。

| 版本 | 输出目录 |
| --- | --- |
| Y 形 | `outputs/y_channel/D1p6_PDMSmin9_ABScore_PLAfixture_bottomless_dovetail_raised3mm_magnets_3dm_P1S04_unsliced/` |
| 三端口直管 | `outputs/straight_channel/D1p6_PDMSmin9_ABScore_PLAfixture_bottomless_dovetail_raised3mm_magnets_3dm_P1S04_internal_end_unsliced/` |
| 双上端口直管 | `outputs/straight_channel/D1p6_PDMSmin9_ABScore_PLAfixture_two_top_ports_raised3mm_magnets_3dm_P1S04_unsliced/` |

每个目录包含两份 `*_project.3mf`、装配/排版 3dm、`assembly/`、`print_parts/`、几何检查图和 `manual_slicing_summary.json`。完整 ZIP 位于 `deliverables/`，文件名包含 `ABScore_PLAfixture` 并以 `_unsliced.zip` 结尾。包内不含已切片 3MF、G-code、临时打印支撑刀路或历史切片报告。

`reference/magnet_slot.3dm` 是用户磁铁槽原件的副本，作为设计来源保留，不要作为额外打印件导入。

## 自行切片

1. 在 Bambu Studio 使用“打开项目”打开 `*_core_ABS_project.3mf`，确认 P1S、0.4 mm 喷嘴和 ABS。血管芯仍采用既有打印姿态；支撑和两层 raft 仅是可编辑的起始设置，尚未生成支撑刀路。
2. 打开 `*_fixture_PLA_project.3mf`，确认 PLA。它使用实际安装的 `Bambu PLA Basic @BBL P1S 0.4 nozzle` 预设，温度、风扇等来自 PLA 预设，没有沿用 ABS 温度。框架和锁键的默认支撑关闭，燕尾支架的默认支撑开启，可按预览修改。
3. 分别检查模型、支撑、孔槽与首层，再由您执行切片。两份工程都可以调整层高、壁数、填充及支撑；当前没有自动生成的打印时间或耗材预测。

浇注时按原燕尾与锁键流程装配，磁铁从上方放入外置槽，配合外部磁性承托板或底面下的对置磁铁固定。无底盒仍须与外部平整底面贴合并处理周边密封。本轮未进行 PLA 配合、磁力、耐温或实际灌注验收。

## 再次生成

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh generate_y_channel_casting.py
./run_python.sh generate_straight_channel_casting.py
./run_python.sh generate_straight_channel_casting.py \
  --config config/straight_channel_two_top_ports_1p6_p1s_abs.yaml
```

文件名中的旧 `p1s_abs` 配置名保留作入口兼容，当前配置正文明确为 ABS 芯、PLA 夹具。默认 `--stage all` 只生成模型、未切片工程和 ZIP。`--stage geometry` 只生成模型；`--stage project` 从已检查的模型导出未切片工程和包。`--stage slice` 已从这两个入口移除，Python 调度接口也拒绝该阶段。

未切片工程只读取本机 Bambu 预设，不需要发现或启动切片引擎。预设根目录可通过 `VASCULAR_BAMBU_PROFILE_ROOT` 指定；服务器使用 `run_server_python.sh` 时沿用服务器本地预设路径。实现集中在 `vascular_processing/channel_manual_delivery.py`，日志和交付检查在 `logs/manual_slice_ABS_PLA_20261012/`。

上一轮全 ABS 已切片交付被移入 `history/before_manual_slice_ABS_PLA_20261012/` 下对应的 `outputs/` 和 `deliverables/`，用于追溯，不属于当前交付。
