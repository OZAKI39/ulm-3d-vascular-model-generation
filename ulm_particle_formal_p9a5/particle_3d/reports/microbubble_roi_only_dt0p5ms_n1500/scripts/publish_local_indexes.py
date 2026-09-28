"""Publish local current-result links only after independent verification passes."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parents[1]
PROJECT=HERE.parents[2]
WORKSPACE=PROJECT.parent
PREVIOUS='microbubble_best_balance_dt0p5ms_n1500'
CURRENT=HERE.name
REL=f'particle_3d/reports/{CURRENT}'
CASE='mean-2p0-mmps-A-ROI-only-balanced-pressure-v1'
SHA='fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12'


def main():
    assert PROJECT==Path('/home/lzy/projects/ulm_particle_formal_p9a5')
    verified=json.loads((HERE/'LOCAL_VERIFICATION.json').read_text())
    result=json.loads((HERE/'data/final_summary.json').read_text())
    assert verified['PASS'] and verified['tracks']==1500
    assert result['flow_sha256']==SHA and result['numerical_gate']
    assert verified['English_animation_and_130percent_projection_verified']
    snapshots=HERE/'audit/before_current_index_update'
    old_workflow=json.loads((snapshots/'CURRENT_WORKFLOW.json').read_text())
    old_results=(snapshots/'CURRENT_RESULTS.md').read_text()
    old_results=old_results.replace('## 最新微泡：best-feasible-balance-v1','## 上一批微泡：best-feasible-balance-v1',1)
    outcomes='、'.join(f'{k}={v}' for k,v in result['outlets'].items())
    statuses='、'.join(f'{k}={v}' for k,v in result['statuses'].items())
    intro=f'''## 最新微泡：ROI-only-balanced-pressure-v1，dt = 0.5 ms

2026-09-28 已完成 1500 条：{outcomes}。终态：{statuses}；执行失败 {result['failed_execution']}。
保持上一批原有粒径分布及固定来源顺序，未按出口补样。数值检查、CUDA 复核、本地独立完整性核验均通过。

- [英文动画（主体放大 30%）、图件与数据入口]({REL}/OPEN_RESULTS.html)
- [中文结果报告]({REL}/MICROBUBBLE_RESULTS_ZH.md)
- [0.5 ms 位置数据]({REL}/data/trajectories_dt0p5ms.npz) · [逐轨迹 CSV]({REL}/data/trajectory_catalog.csv) · [出口统计]({REL}/data/outlet_summary.csv)
- [1500 条原始轨迹及审计]({REL}/tracks/)
- [本地独立核验]({REL}/LOCAL_VERIFICATION.json) · [复现说明]({REL}/REPRODUCE.md)

8 CPU 进程执行有限尺寸积分，RTX 4090 用于 CUDA FP64 场/轨迹复核和 EGL 动画渲染；视频编码为 CPU libx264。
动画按轨迹年龄对齐，不代表同时注入；1.30 倍为相对原始显示比例，未重复叠加。所有既有结果继续保留。

'''
    (PROJECT/'CURRENT_RESULTS.md').write_text(old_results.replace('# 当前结果入口\n\n','# 当前结果入口\n\n'+intro,1))
    readme=(snapshots/'README.md').read_text().replace(PREVIOUS,CURRENT).replace('best-feasible-balance-v1','ROI-only-balanced-pressure-v1')
    (PROJECT/'README.md').write_text(readme)
    workflow=dict(old_workflow)
    workflow.update(schema='CURRENT_WORKFLOW_20260928_ROI_ONLY_MB1500',particle_report=REL,
        flow_sha256=SHA,flow_loader=REL+'/scripts/campaign.py:setup',
        source_cohort=REL+'/data/cohort.json',formal_contract=REL+'/data/contract.json',
        verification=REL+'/LOCAL_VERIFICATION.json',
        server_report='/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928',
        preserved_previous_batch=old_workflow)
    (PROJECT/'CURRENT_WORKFLOW.json').write_text(json.dumps(workflow,indent=2,ensure_ascii=False)+'\n')
    agents=(snapshots/'AGENTS.md').read_text().replace(PREVIOUS,CURRENT).replace('best-feasible-balance-v1 flow','ROI-only-balanced-pressure-v1 flow')
    (PROJECT/'AGENTS.md').write_text(agents)
    active=(snapshots/'ACTIVE_VASCULAR_WORKFLOW.md').read_text()
    active=active.replace('mean-2p0-mmps-A-best-feasible-balance-v1',CASE).replace('当前 FEM：best-feasible-balance-v1','当前 FEM：ROI-only-balanced-pressure-v1')
    active=active.replace('fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2',SHA).replace('steady_flow_mean_2p0_mmps_A_best_feasible_balance.vtu','steady_flow_mean_2p0_mmps_A_ROI_only_balanced_pressure.vtu')
    (WORKSPACE/'ACTIVE_VASCULAR_WORKFLOW.md').write_text(active)
    server=(snapshots/'CURRENT_SERVER_PATHS.md').read_text()
    server=server.replace('| 当前 FEM：best-feasible-balance-v1 |','| 上一版 FEM：best-feasible-balance-v1 |').replace('| 当前微泡：1500 条，dt=0.5 ms |','| 上一批微泡：1500 条，dt=0.5 ms |')
    marker='|---|---|\n'
    new_rows=f'| 当前 FEM：ROI-only-balanced-pressure-v1 | `/workspace/flow_roi_only_balance_20260928/{CASE}` |\n| 当前微泡：1500 条，dt=0.5 ms | `/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928` |\n'
    server=server.replace(marker,marker+new_rows,1)
    (WORKSPACE/'CURRENT_SERVER_PATHS.md').write_text(server)
    (HERE/'RUNNING_STATUS.md').write_text(f'# 已完成并通过本地独立核验\n\n1500/1500，dt=0.5 ms；{outcomes}。终态：{statuses}。\n\n入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)。运动积分、CUDA 复核及英文 1.30 倍动画均已完成。\n')
    print('Published verified ROI-only microbubble results and preserved prior result links')


if __name__=='__main__':main()
