"""Create a read-only manual data retention review; never move or delete data."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Old names do not establish obsolescence. Some reference files remain necessary.
GROUPS = {
    'current': ('继续保留', '最新 ROI-only 成果', '最新流场、1500 条微泡、配置、计算日志和图表。'),
    'geometry': ('继续保留', '最新方案的原始血管与几何', '后续继续开发所需的原始样本、ROI 和表面；不是旧方法结果。'),
    'network': ('继续保留；可另做细分', '最新方案仍使用的网络输入', '文件名带 H0/v1/v2，但 network.py 直接读取其中的图、端口、参考状态及压力设计证据；不能整包删除。'),
    'reference': ('继续保留', '最新程序共用的参考输入', '包含原网格、原参考场和粒径/来源身份；flow_field.load_case 会校验/读取旧参考场，不能因名称旧就删除。'),
    'alternative': ('待你确认', '另一套压力方案 best-feasible', '与最新方案并非同一结果；属于对照资料，不认定为无效。删除会使相应测试、比较图和全量核验失效。'),
    'h0': ('待你确认；先保护一个小文件', '早期 H0 正式结果', '新 ROI-only 小测试仍读取 results/network_h0/flow/run/solver.xml；整包删除前须把该文件独立保留并更新测试引用。'),
    'baseline': ('待你确认', '原始 2 mm/s 算例结果', '独立算例输出可审查；inputs/reference/fem 中同名参考输入仍须保留，两者不能一起删。'),
    'stage_q': ('待你确认；先保留环境记录', '早期 Stage-Q 结果', '完整旧场可审查，但 run/runtime_manifest.json 与 gpu_environment_manifest.json 是原求解软件环境依据，建议保留小型记录。'),
    'validation': ('保留验证依据；大文件可单列审查', '现方法的已有验证资料', '圆管、步长、网格和壁面量记录仍支持当前方法，不能仅按日期认定过时；删大结果前保留验证报告和文件身份。'),
    'history': ('待你确认', '精选历史说明和旧执行脚本', '不参与当前计算；体积很小，建议留作追溯，也可手动决定删除。'),
    'other': ('继续保留', '配置、源码或安装依据', '原迁移清单中的其他必要材料；不是本轮拟删的旧结果。'),
}


def group(path):
    for prefix, category in (
        ('results/roi_only/', 'current'), ('inputs/geometry/', 'geometry'),
        ('inputs/network/', 'network'), ('inputs/reference/', 'reference'),
        ('results/best_feasible/', 'alternative'), ('results/network_h0/', 'h0'),
        ('results/mean_2mmps/', 'baseline'), ('results/stage_q/', 'stage_q'),
        ('validation/', 'validation'), ('archive/', 'history'),
    ):
        if path.startswith(prefix):
            return category
    return 'other'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows = json.loads((ROOT / 'docs/audit/data_manifest.json').read_text())
    summaries = defaultdict(lambda: {'files': 0, 'bytes': 0, 'top_paths': set()})
    with (output / 'FILES.csv').open('x', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['分类', '建议', '新工程内路径', '原位置', '字节数', '原核验SHA256', '你的决定'])
        for row in rows:
            category = group(row['path'])
            decision, title, reason = GROUPS[category]
            writer.writerow([title, decision, row['path'], row['source'], row['size'], row['sha256'], ''])
            summaries[category]['files'] += 1
            summaries[category]['bytes'] += row['size']
            summaries[category]['top_paths'].add('/'.join(row['path'].split('/')[:2]))
    items = []
    for key, (decision, title, reason) in GROUPS.items():
        summary = summaries[key]
        items.append({'id': key, 'title': title, 'recommendation': decision, 'reason': reason,
                      'files': summary['files'], 'bytes': summary['bytes'],
                      'top_paths': sorted(summary['top_paths']), 'user_decision': ''})
    with (output / 'GROUPS.csv').open('x', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['资料组', '文件数', '大小MB', '建议', '注意事项', '你的决定'])
        for item in items:
            writer.writerow([item['title'], item['files'], round(item['bytes'] / 1e6, 2), item['recommendation'], item['reason'], ''])
    report = {'current_workflow': 'roi_only', 'scope': 'Previously verified migrated data; excludes new venv/caches and old workspaces outside this project.',
              'deleted': [], 'moved': [], 'groups': items,
              'warning': 'This is a review inventory, not a deletion script. After an approved deletion, indexes, comparison links, tests and active integrity manifests must be updated without rewriting original audit evidence.'}
    (output / 'inventory.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    lines = ['# 最新资料保留与旧资料人工审查', '',
             '按当前工程的默认 ROI-only 方案作为最新工作流。清单已完成分类，文件未移动、未删除。你可在 GROUPS.csv 的“你的决定”列填“保留”或“同意删除”。', '',
             'best-feasible 是另一方案，不等于无效旧结果。带旧名称的共享输入仍可能是最新程序必需的。', '',
             '| 资料组 | 文件数 | 大小 MB | 建议 |', '|---|---:|---:|---|']
    for item in items:
        lines.append(f"| {item['title']} | {item['files']} | {item['bytes']/1e6:.2f} | {item['recommendation']} |")
    lines += ['', '这里的大小是新工程中已核验副本的实际文件字节总数，不含原目录和 Python 环境。删除新副本不会自动释放旧目录所占空间。', '', '## 删除前必须看', '']
    for item in items:
        lines.append('- **' + item['title'] + '**：' + item['reason'])
    lines += ['', '当前程序还有面向全部旧算例的检查和比较入口；手动删完旧结果会让它们找不到文件。确认删除范围后，应一起更新当前索引和这些入口，原始迁移审计仍保留原文。', '',
              '最新小测试本身仍用到 H0 的 solver.xml；不要直接删除整个 H0 文件夹。旧名称输入中的 mesh、流场参考、网络端口、SonoVue 数据和源码身份证据均标为保留。', '',
              '逐文件路径、原位置、大小和此前已核验的 SHA 在 [FILES.csv](FILES.csv)；便于人工决定的汇总在 [GROUPS.csv](GROUPS.csv)。', '',
              '旧项目外部目录的清理仍使用 [原目录清单](../CLEANUP_REVIEW_ZH.md)。本清单不授权删除共享目录、兼容链接或其他项目。', '']
    (output / 'README.md').write_text('\n'.join(lines))
    print(json.dumps({'files': len(rows), 'bytes': sum(r['size'] for r in rows), 'output': str(output), 'deleted': 0}, ensure_ascii=False))


if __name__ == '__main__':
    main()
