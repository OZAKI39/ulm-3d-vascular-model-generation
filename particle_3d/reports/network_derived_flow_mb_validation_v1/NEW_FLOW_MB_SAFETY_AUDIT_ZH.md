# NEW 30 泡安全审核

最终安全门槛 PASS。穿透、采样态 handoff 违规、连续证书违规、NaN/Inf、未分类 solver corruption、solver fail 均为 0；无 inlet escape。所有 30 个结果都有终态记录，没有删除不完成轨迹。

共检查 12147 个顶层连续证书，递归得到 13025 个原始 support-plane-minus-h_lower 叶证书。每个 union 的子区间无重叠空隙地覆盖 [0,1]，且声明 held-velocity path 未改；每个叶证书的 minimum_g_nf_bound 均满足原 roundoff budget。采样态 gap≥−roundoff 与 g_nf≥−roundoff 分别复核；roundoff=2.009411463123824e-17 m。

## 唯一 stationary：ID 15

OLD/NEW 中均为 D=2.489958478 µm。最终 32 个已接受求解具有 rank=3 的三个独立接触约束，乘子非负，平移速度小于原 KKT velocity budget。末态坐标约 (92.886303,46.631184,114.771900) µm，gap≈2 nm；NEW 局部流速 6.081224 mm/s。故支持当前刚性有限尺寸接触静止，不能声称流体停滞、软件失败、体内捕获证据或真实微泡必然无法通过。无半径替代实验。

## 确定性与科学保护

1/3 workers 的 ID 3、15、22 全 samples 与解压事件日志逐位一致；永久测试又与 6 workers 正式结果逐位比较。12 项测试全部通过。首次测试的证书 schema KeyError 是新观察代码对 union 证书支持不足，原始失败日志已保留于 `logs/initial_tests.log/xml`，修复限于本轮新测试与本轮新审计器。

本地包内原有 206 个文件、服务器源码快照 114 个文件均 SHA 未改变。NEW/OLD formal VTU、复用 OLD 轨迹及原报告均未覆盖，新增报告、adapter 和测试位于独立新目录。未执行任何科学源码变更、FEM 修正或新增 500 泡任务。

## 证据

[所有末态及原求解器诊断](data/saved_state_safety_audit.json)、[NEW 指标](data/new_paired_metrics.json)、[12 tests](logs/tests.log)、[science protection](data/particle_science_protection_final.json)、[remote snapshot](logs/remote_final_snapshot_verification.json)、[worker determinism](data/worker_count_determinism_subset.json)。

NEW P1 局部截面误差 max=3.454212%、RMS=1.789924% 仍是本实验已知模型限制。建议 `READY_FOR_LARGE_SAMPLE_REVIEW`，需人工审核后再决定大样本。
