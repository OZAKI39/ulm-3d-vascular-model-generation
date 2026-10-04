"""Build the local reading page from the reviewed Chinese report."""
from pathlib import Path
from markdown_it import MarkdownIt

root = Path(__file__).resolve().parent
body = MarkdownIt('commonmark').enable('table').render((root / 'CAUSE_ANALYSIS_ZH.md').read_text())
page = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>BraVa 微泡失败与速度急变分析</title>
<style>
body{max-width:1040px;margin:40px auto;padding:0 22px;font:17px/1.8 system-ui,sans-serif;color:#18232d;background:#fafbfc}
h1{font-size:30px;line-height:1.4}h2{font-size:23px;margin-top:2em}
img{max-width:100%;height:auto;background:white}a{color:#176597}
table{border-collapse:collapse;width:100%;font-size:15px}td,th{padding:8px;border:1px solid #ccd3db;text-align:left}
pre{overflow:auto;padding:16px;background:#eaf0f5}code{font-size:.88em;overflow-wrap:anywhere}
</style><body>'''+body+'</body></html>\n'
(root / 'OPEN_RESULTS.html').write_text(page)
print(root / 'OPEN_RESULTS.html')
