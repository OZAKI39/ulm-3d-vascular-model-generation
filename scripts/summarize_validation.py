"""Summarize actual local receipts without running any scientific calculation."""
import datetime
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from human_vessel_flow.paths import ROOT, read_json, write_json

def receipt(path, label, scope, status=None):
    p=ROOT/path
    if not p.exists():return dict(label=label,status='NOT_RUN',scope=scope,evidence=path)
    data=read_json(p)
    value=status or ('PASS' if str(data.get('status','')).startswith('PASS') else 'FAIL')
    return dict(label=label,status=value,scope=scope,evidence=path)

def main():
    rows=[
        receipt('validation/FINAL_IMMUTABLE_DATA.json','冻结资料逐字节身份','两套场、全部旧群体、输入、历史证据与原生源码'),
        receipt('docs/provenance/ORIGINAL_FILES_UNCHANGED.json','原来源文件未改动','3439 个原路径的 SHA 复查'),
        receipt('validation/ORIGINAL_WORKTREE_STATUS.json','原 Git 工作树','相同参数核对 HEAD、状态数量和既有源码变更列表'),
        receipt('validation/GEOMETRY_AND_POSES.json','几何、单位、两姿态','包含完整源表面坐标往返与速度大小保持'),
        receipt('validation/read_only_v1/mesh.json','原网格检查','节点、单元、边界、质量和封闭性'),
        receipt('validation/read_only_v1/cfd_wss.json','冻结 CFD / WSS / checkpoint','只读复核；没有新求解，不认证近壁精度'),
        receipt('validation/read_only_v1/population.json','旧群体状态和账本','保留 141 观察截止、2 数值失败及 false 科学完整性状态'),
        receipt('validation/reference_reproduction/REVIEW.json','V2 已存误差分解复算','原科学结论仍为 MIXED_FIELD_ERROR'),
        receipt('validation/INDEPENDENCE.json','源码与数据独立性','隔离 Python，拒绝旧工作区读取，实际导入与数据载入'),
        receipt('validation/real_gpu_field/RESULT.json','真实 BraVa 场 CPU/GPU 求值','2146 位置；无轨迹积分，不能证明局部物理准确'),
        receipt('validation/MESHER_ENVIRONMENT.json','网格软件加载','SimVascular 内嵌 Python 导入 TetGen；没有生成新网格'),
        receipt('validation/NATIVE_ENVIRONMENT.json','原生求解器构建和加载','实际编译、动态库和无输入启动；没有执行 CFD',
                'PASS' if (ROOT/'validation/NATIVE_ENVIRONMENT.json').exists() and all(x['status']=='PASS' for x in read_json(ROOT/'validation/NATIVE_ENVIRONMENT.json')['checks'].values()) else 'FAIL'),
        receipt('validation/CURRENT_LINKS.json','当前报告和可视化链接','新报告、新索引、新 HTML；历史 HTML 的旧资源链接另列'),
    ]
    junit=ROOT/'validation/small_tests/junit.xml';test_stats={}
    if junit.exists():
        suites=ET.parse(junit).getroot();cases=suites.findall('.//testcase')
        failures=sum(c.find('failure') is not None for c in cases);errors=sum(c.find('error') is not None for c in cases)
        skipped=sum(c.find('skipped') is not None for c in cases)
        test_stats=dict(total=len(cases),passed=len(cases)-failures-errors-skipped,failed=failures,errors=errors,skipped=skipped)
        rows.append(dict(label='原 CPU/GPU 小型测试及工程保护测试',status='FAIL' if failures or errors else 'PASS',
                         scope=str(test_stats)+'；人工小场景，没有正式 BraVa 群体',evidence='validation/small_tests/pytest.log'))
    else:rows.append(dict(label='原 CPU/GPU 小型测试及工程保护测试',status='NOT_RUN',scope='等待实际 JUnit 记录',evidence='validation/small_tests/pytest.log'))
    rows.extend([
        dict(label='新正式 BraVa CFD / 新网格',status='NOT_RUN',scope='本轮禁止启动；仅准备了计划',evidence='validation/prepared_case_example/PLAN.json'),
        dict(label='新固定时间窗正式群体',status='NOT_RUN',scope='协议与小测试已另列；没有新正式批次',evidence='validation/prepared_population_example/PLAN.json'),
        dict(label='V2 局部流场精度',status='FAIL',scope='保留原科学判定 REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC',evidence='diagnostics/reference_pipe/FIELD_VALIDATION.json'),
        dict(label='P1 轨迹敏感性 / 收敛认证',status='FAIL',scope='保留 P1_FIELD_TRAJECTORY_IMPACT_LARGE；并非本轮新增失败',evidence='diagnostics/reference_pipe/P1_FIELD_PARTICLE_SENSITIVITY.json'),
        dict(label='P2/P1 正式可用性',status='BLOCKED',scope='原资源限制停止，且曲面压力检查仍有失败；本轮未恢复',evidence='diagnostics/p2p1/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md'),
        dict(label='V3 新网格和新场',status='NOT_RUN',scope='保留提案身份',evidence='diagnostics/reference_pipe/PROPOSED_V3_MESH_POLICY.json'),
        dict(label='可选旧 LAMMPS 接口',status='NOT_RUN',scope='不是当前 BraVa 群体所用后端；未安装测试',evidence='docs/NATIVE_ENVIRONMENT_ZH.md'),
    ])
    result=dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                project_refactor_scope='Local independent engineering refactor; frozen scientific states preserved',
                checks=rows,test_counts=test_stats,new_formal_CFD_runs=0,new_formal_population_runs=0,
                original_files_deleted=0,remote_connections=0,github_pushes=0,
                overall_scientific_accuracy='NOT_CERTIFIED')
    write_json(ROOT/'validation/DELIVERY_STATUS.json',result)
    md=['# 本轮实际验证结果','',
        '以下 PASS 指对应的具体检查，不是整套流场物理精度通过。FAIL/BLOCKED 的旧科学结论完整保留。','',
        '| 检查 | 状态 | 实际范围 / 证据 |','|---|---|---|']
    for row in rows:
        md.append('| '+row['label']+' | **'+row['status']+'** | '+row['scope']+'；[记录](../'+row['evidence']+') |')
    md.extend(['','完整测试使用本地实际 CUDA 设备；是否有跳过以 JUnit 数字为准。初次缺少依赖的收集失败、CPU 先行测试中的 GPU 跳过、安装超时及重试都保存在 `validation/logs/`，最终结果使用完整环境下的独立测试目录。',
               '', '本轮没有完整 CFD 数值回归；原生二进制能启动，只是环境与接口检查。也没有把旧 1565 条结果当成新协议的大批量运行。',
               '', '实际命令入口见 [README](../README.md)，新输出目录已经存在时会拒绝覆盖；复查请选择新目录。'])
    (ROOT/'docs/VALIDATION_SUMMARY.md').write_text('\n'.join(md)+'\n')
    print(json.dumps(dict(tests=test_stats,checks=len(rows)),indent=2))

if __name__=='__main__':main()
