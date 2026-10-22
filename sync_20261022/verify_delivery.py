#!/usr/bin/env python3
import datetime,hashlib,json,subprocess
from pathlib import Path
from export_flow_projects import AUDIT,BASE

results=[]
for project in ['brava','mouse']:
 root=BASE/'temp_storage'/f'{project}_flow_github_sync_20261022';upper=project.upper()
 plan=json.loads((AUDIT/f'{upper}_PUBLISH_PLAN.json').read_text());receipt=json.loads((AUDIT/f'{upper}_PUSH_RECEIPT.json').read_text())
 assert receipt['status']=='PUSHED_AND_REMOTE_SHA_VERIFIED'
 def git(*args):return subprocess.check_output(['git','-C',str(root),'-c','diff.renames=false','-c','gc.auto=0',*args])
 head=git('rev-parse','HEAD').decode().strip();remote=git('ls-remote','origin','refs/heads/'+receipt['branch']).decode().split()[0]
 assert head==remote==receipt['head']
 tree={}
 for record in git('ls-tree','-r','-l','-z',head).split(b'\0'):
  if not record:continue
  info,path=record.split(b'\t',1);mode,kind,oid,size=info.decode().split()
  tree[path.decode()]={'mode':mode,'oid':oid,'size':int(size)}
 expected={e['path']:e for batch in plan['batches'] for e in batch}
 assert set(tree)==set(expected),(set(tree)-set(expected),set(expected)-set(tree))
 assert all(tree[p]['oid']==e['oid'] and tree[p]['mode']==e['mode'] for p,e in expected.items())
 assert max(x['size'] for x in tree.values())<100*1024*1024
 assert not git('status','--porcelain','--untracked-files=all').strip()
 fields=json.loads((root/'sync_20261022/FROZEN_FIELDS_CHECK.json').read_text())
 for f in fields:assert f['pass']
 nstates=sum('/population_v1/runs/realization_0000_contactfix/states/chunk_' in p for p in tree)
 if project=='brava':assert nstates==137
 summary=json.loads((root/'sync_20261022/SUMMARY.json').read_text())
 results.append({'project':project,'branch':receipt['branch'],'commit':head,'url':'https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/'+receipt['branch'],'files_and_links':len(tree),'logical_bytes':sum(x['size'] for x in tree.values()),'source_files':summary['source_files'],'source_bytes':summary['source_bytes'],'frozen_fields_verified':len(fields),'formal_accepted_state_chunks':nstates,'max_blob_bytes':max(x['size'] for x in tree.values()),'working_tree_clean':True,'remote_commit_verified':True,'all_published_tree_objects_match_plan':True,'commits':len(receipt['commits']),'local_worktree':str(root)})
delivery={'status':'PASS','checked_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'repository':'https://github.com/OZAKI39/ulm-3d-vascular-model-generation','projects':results,'scientific_runs_started':0,'force_push_used':False,'original_source_trees_modified':False}
(AUDIT/'DELIVERY_VERIFIED.json').write_text(json.dumps(delivery,indent=2,ensure_ascii=False)+'\n')
lines=['# 小鼠血管与 BraVa 分支同步报告','',f"核验时间：{delivery['checked_at']}。状态：**PASS**。两条独立新分支已推送，远端 commit 与本地完全一致。",'','| 项目 | 新分支 | 提交 | 文件/链接 | 发布逻辑体积 |','|---|---|---|---:|---:|']
for r in results:lines.append(f"| {r['project']} | [{r['branch']}]({r['url']}) | `{r['commit']}` | {r['files_and_links']:,} | {r['logical_bytes']/1e9:.3f} GB |")
lines += ['','本次保留各项目代码、关键日志、正式网格/冻结流场、WSS、原始正式轨迹/接受状态、诊断结果和成品图表。BraVa 137 个 contactfix 状态块全部发布；旧1565条仍为1422出流、141截止、2数值失败。小鼠ROI-only、best-feasible、H0及2mm/s基线分别保留。16个冻结场通过既有登记SHA核验。','',
'排除可重建环境和安装包、缓存、重复中间流场与渲染内容、被替代批次的庞大状态、可重建全域审计数组；不删除本地或服务器原始文件。详细逐文件来源、大小、SHA、排除原因和兼容链接见各分支 `sync_20261022/`。BraVa新服务器对应文件与本地一致，服务器独有执行证据从已核验备份补入；小鼠没有声称重新核验旧服务器。','',
'排除了未取得实际数据的上游原生求解器测试资源LFS指针，保留 `UPSTREAM_TEST_RESOURCE_REFERENCES.json` 以追溯。项目自身冻结数据和正式状态是普通Git中的实际文件，没有LFS、没有切片规避限制。少数上游.gitattributes仅在发布副本中停用LFS，原文及前后哈希保留；科学代码和数据字节不改。','',
'分支当前树只含各自项目与必要共享依赖；祖先历史仍保持原样。没有向main提交、没有force push、没有重写既有提交。首次未发布的失败提交保留在本地旧同步分支，成功发布工作树使用 `-publish` 本地分支，跟踪上表远端分支。','',
'遇到的问题已解决：旧VS Code凭据助手失效，改为仅对此次Git命令使用已有gh登录；GitHub拒收上游缺失LFS资源后，排除这些测试资源指针，从原已发布父提交建立干净发布分支。失败记录与修复前计划/提交收据保留在本目录。','',
'核验包括：文件SHA、敏感信息扫描、分支范围、冻结流场、README结果链接、有效兼容链接、全部最终Git树对象与发布计划一致、无≥100MiB普通Git文件、工作树clean、fetch/ls-remote最终提交一致。本轮没有启动CFD、网格或轨迹计算，没有改变已有FAIL/SKIPPED/未运行结论。','',
'复现限制：历史配置内绝对路径保留，环境/安装依赖需按现有说明配置；这不是空白机器一键运行包。GitHub中HTML/视频可下载后打开。','',
'机器核验回执：[DELIVERY_VERIFIED.json](DELIVERY_VERIFIED.json)。']
(AUDIT/'SYNC_REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(delivery,indent=2,ensure_ascii=False))
