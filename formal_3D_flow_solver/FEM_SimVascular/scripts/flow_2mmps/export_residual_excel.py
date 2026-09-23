"""Export the current residual panels as separate, editable Excel workbooks.

Reads audited CSV data only; does not run the solver or change existing figures.
Run with the project Python environment (requires XlsxWriter).
"""
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile, ZIP_DEFLATED

import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "flow_cases/mean-2p0-mmps/field_diagnostics"
DEST = OUT / "excel"
BLUE, ORANGE, GRAY = "#19496b", "#b97422", "#6c757e"
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "c": "http://schemas.openxmlformats.org/drawingml/2006/chart"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def typed(value):
    if value == "":
        return None
    if value in ("True", "False"):
        return value == "True"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    try:
        return float(value)
    except ValueError:
        return value


def read_csv(name):
    with (OUT / "data" / name).open(newline="") as stream:
        return [{k: typed(v) for k, v in row.items()} for row in csv.DictReader(stream)]


class Book:
    def __init__(self, filename, title, description, source_hashes):
        self.path = DEST / filename
        self.title = title
        self.description = description
        self.sources = source_hashes
        self.book = xlsxwriter.Workbook(self.path)
        self.book.set_properties({"title": title, "subject": description,
                                  "comments": "Audited residual data; editable native Excel chart."})
        self.head = self.book.add_format({"bold": True, "font_color": "white", "bg_color": BLUE,
                                         "text_wrap": True, "valign": "vcenter", "border": 0})
        self.scientific = self.book.add_format({"num_format": "0.000000000000E+00"})
        self.decimal = self.book.add_format({"num_format": "0.00000000"})
        self.integer = self.book.add_format({"num_format": "0"})
        self.wrap = self.book.add_format({"text_wrap": True, "valign": "top"})
        self.expected = {}
        self.charts = []
        self.plot_rows = 0

    def table(self, name, headers, rows, definitions=None, decimal_columns=()):
        sheet = self.book.add_worksheet(name)
        sheet.freeze_panes(1, 1)
        sheet.hide_gridlines(2)
        sheet.set_row(0, 44)
        sheet.set_column(0, len(headers) - 1, 26)
        sheet.write_row(0, 0, headers, self.head)
        sheet.autofilter(0, 0, len(rows), len(headers) - 1)
        for j, header in enumerate(headers):
            if definitions and header in definitions:
                sheet.write_comment(0, j, definitions[header])
        for i, row in enumerate(rows, 1):
            for j, value in enumerate(row):
                if value is None:
                    continue  # Missing groups/tails remain blank, never zero.
                if isinstance(value, bool):
                    sheet.write_boolean(i, j, value)
                elif isinstance(value, int):
                    sheet.write_number(i, j, value, self.integer)
                elif isinstance(value, float):
                    assert math.isfinite(value)
                    sheet.write_number(i, j, value,
                                       self.decimal if j in decimal_columns else self.scientific)
                else:
                    sheet.write_string(i, j, value)
        self.expected[name] = [headers] + rows
        if name == "Plot_Data":
            self.plot_rows = len(rows)
            sheet.set_tab_color(BLUE)
        return sheet

    def raw(self, name, records):
        headers = list(records[0])
        self.table(name, headers, [[r[h] for h in headers] for r in records])

    def chart(self, xlabel, ylabel, *, kind="scatter", logarithmic=True,
              ymin=None, ymax=None, xmax=None):
        chart = self.book.add_chart({"type": kind, **({"subtype": "straight"} if kind == "scatter" else {})})
        font = {"name": "Arial", "size": 11}
        chart.set_title({"name": self.title, "name_font": {"name": "Arial", "size": 15}})
        chart.set_chartarea({"fill": {"color": "white"}, "border": {"none": True}})
        chart.set_plotarea({"fill": {"color": "white"}, "border": {"none": True}})
        chart.set_legend({"position": "bottom", "font": font})
        xaxis = {"name": xlabel, "name_font": font, "num_font": font, "num_format": "0",
                 "major_tick_mark": "outside", "line": {"color": GRAY},
                 "major_gridlines": {"visible": False}}
        if kind == "scatter":
            xaxis.update({"min": 0, "max": xmax, "major_unit": 100 if xmax and xmax > 200 else 10})
        else:
            xaxis["interval_unit"] = 10
        yaxis = {"name": ylabel, "name_font": font, "num_font": font,
                 "major_gridlines": {"visible": True, "line": {"color": "#e3e7eb", "width": 0.5}},
                 "major_tick_mark": "outside", "line": {"color": GRAY},
                 "num_format": "0E+00" if logarithmic else "0"}
        if logarithmic:
            yaxis["log_base"] = 10
        if ymin is not None:
            yaxis["min"] = ymin
        if ymax is not None:
            yaxis["max"] = ymax
        chart.set_x_axis(xaxis)
        chart.set_y_axis(yaxis)
        chart.show_blanks_as("gap")
        chart.set_size({"width": 1150, "height": 660})
        sheet = self.book.add_worksheet("Example_Chart")
        sheet.hide_gridlines(2)
        sheet.insert_chart("A1", chart)
        sheet.set_column("A:A", 155)
        sheet.write("A36", "可编辑示例图：数据与原图相同；论文级注释和排版请按需调整。详见 Notes。", self.wrap)
        sheet.set_tab_color(ORANGE)
        self.charts.append((chart, logarithmic))
        return chart

    def series(self, chart, column, label, color, *, sheet="Plot_Data", first=1, last=None,
               marker=None, connect=True, dashed=False):
        if last is None:
            last = self.plot_rows
        line = {"color": color, "width": 1.3} if connect else {"none": True}
        if dashed:
            line["dash_type"] = "dash"
        chart.add_series({"name": label, "categories": [sheet, first, 0, last, 0],
                          "values": [sheet, first, column, last, column], "line": line,
                          "marker": marker or {"type": "none"}})

    def notes(self, specific, definitions):
        rows = [("对应图", self.description),
                ("使用方式", "Plot_Data 第一行是列名，第二行起为数值。按下列 X/Y 列插入图表；Example_Chart 可直接编辑、复制。"),
                ("数字与空值", "残差以数值存储，不是科学计数法文本。保留源 CSV 有效精度（Excel 约 15 位有效数字）；空白代表无该组数据，不能填 0。"),
                ("物理意义", "归一化残差和残差比值无量纲；残差大小不能直接解释为 CFD 解的物理误差。未证明存在浮点精度极限。")]
        rows.extend(specific)
        rows.extend((f"列：{key}", value) for key, value in definitions.items())
        rows.extend((f"来源 SHA256：{key}", value) for key, value in self.sources.items())
        sheet = self.table("Notes", ["Item", "说明 / Definition"], [list(r) for r in rows])
        sheet.set_column(0, 0, 44, self.wrap)
        sheet.set_column(1, 1, 112, self.wrap)
        for i in range(1, len(rows) + 1):
            sheet.set_row(i, 48)

    def close_and_verify(self):
        self.book.close()
        numeric = 0
        blank = 0
        with ZipFile(self.path) as archive:
            assert archive.testzip() is None
            shared = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(item.itertext()) for item in shared]
            workbook = ET.fromstring(archive.read("xl/workbook.xml"))
            sheet_names = [s.attrib["name"] for s in workbook.find("s:sheets", NS)]
            for name, rows in self.expected.items():
                index = sheet_names.index(name) + 1
                root = ET.fromstring(archive.read(f"xl/worksheets/sheet{index}.xml"))
                cells = {c.attrib["r"]: c for c in root.findall(".//s:sheetData/s:row/s:c", NS)}
                for i, row in enumerate(rows):
                    for j, value in enumerate(row):
                        address = xl_rowcol_to_cell(i, j)
                        cell = cells.get(address)
                        if value is None:
                            assert cell is None or cell.find("s:v", NS) is None
                            blank += 1
                            continue
                        assert cell is not None, (name, address)
                        saved = cell.find("s:v", NS).text
                        if isinstance(value, bool):
                            assert cell.attrib["t"] == "b" and int(saved) == int(value)
                        elif isinstance(value, (int, float)):
                            assert cell.attrib.get("t", "n") == "n"
                            assert math.isclose(float(saved), value, rel_tol=2e-15, abs_tol=0), (name, address)
                            numeric += 1
                        else:
                            assert strings[int(saved)] == value, (name, address)
            chart_files = [n for n in archive.namelist() if re.fullmatch(r"xl/charts/chart\d+.xml", n)]
            assert len(chart_files) == len(self.charts) == 1
            chart = ET.fromstring(archive.read(chart_files[0]))
            log_nodes = chart.findall(".//c:logBase", NS)
            assert bool(log_nodes) == self.charts[0][1]
            assert all(x.attrib["val"] == "10" for x in log_nodes)
            blanks_mode = chart.find(".//c:dispBlanksAs", NS)
            # OOXML defaults to gaps; XlsxWriter omits the default element.
            assert blanks_mode is None or blanks_mode.attrib["val"] == "gap"
            assert "#REF!" not in archive.read(chart_files[0]).decode()
            series = chart.findall(".//c:ser", NS)
            assert len(series) == len(self.charts[0][0].series)
        return {"file": self.path.name, "sha256": sha(self.path), "plot_data_rows": self.plot_rows,
                "numeric_cells_verified": numeric, "blank_cells_verified": blank,
                "editable_chart_series": len(series), "y_axis": "log10" if self.charts[0][1] else "linear",
                "numeric_roundtrip_relative_tolerance": 2e-15, "all_pass": True}


def main():
    DEST.mkdir(exist_ok=True)
    source_names = ["data/nonlinear_step_audit.csv", "data/residual_linear_solves.csv",
                    "data/residual_true_monitor.csv", "NONLINEAR_RESIDUAL_AUDIT.json", "COMPUTE_VALIDATION.json"]
    source_hashes = {name: sha(OUT / name) for name in source_names}
    steps = read_csv("nonlinear_step_audit.csv")
    solves = read_csv("residual_linear_solves.csv")
    selected = [solves[k]["linear_solve_index"] for k in (0, 83, 166)]
    histories = {index: [] for index in selected}
    with (OUT / "data/residual_true_monitor.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            index = int(row["linear_solve_index"])
            if index in histories:
                histories[index].append({k: typed(v) for k, v in row.items()})
    audit = json.loads((OUT / "NONLINEAR_RESIDUAL_AUDIT.json").read_text())
    compute = json.loads((OUT / "COMPUTE_VALIDATION.json").read_text())
    # Find the already audited linear tolerances without hard-coding them.
    linear = next(value for value in compute.values() if isinstance(value, dict) and "rtol" in value and "atol" in value)
    tol, rtol, atol = audit["nonlinear_tolerance"], linear["rtol"], linear["atol"]
    assert [r["step"] for r in steps] == list(range(1, 72))
    assert [r["linear_solve_index"] for r in solves] == list(range(1, 168))
    for index, rows in histories.items():
        assert [r["iteration"] for r in rows] == list(range(solves[index - 1]["linear_iterations"] + 1))
        assert all(r["true_relative_residual"] > 0 and r["finite"] for r in rows)
        assert rows[-1]["true_relative_residual"] == solves[index - 1]["final_true_relative_residual"]
    for row in steps:
        assert math.isclose(row["reduction_orders"], math.log10(row["start_norm"] / row["end_norm"]), rel_tol=1e-14)
        assert math.isclose(row["start_global"], row["start_norm"] / audit["reference_norm"], rel_tol=1e-14)
    nonlinear_notes = [
        ("固定参考值", f"R_ref = 第 1 时间步第 1 次非线性记录的范数 = {audit['reference_norm']:.14g}；全程固定，不逐步重新归一化。"),
        ("Start / Last", "每一步首次 / 最后一次记录的非线性残差（来自该次线性修正的初始 KSP 范数）。Last 不是最后一次更新之后额外测量的残差。"),
        ("非线性停止条件", f"k >= {audit['maximum_iterations']}，或 k >= {audit['minimum_iterations']} 且（R/R_ref <= {tol:g} 或 R/R_step_start <= {tol:g}）。本次没有达到迭代上限。"),
        ("第 4、23 步", "非线性迭代次数分别由 4 减为 3、由 3 减为 2；均满足停止条件。")]
    results = []

    # A: absolute run-wide normalization of the first record in each step.
    book = Book("A_Global_Transient_Evolution.xlsx", "(A) Global transient evolution",
                "residual_convergence (A) / residual_nonlinear 左图", source_hashes)
    definitions = {"Step": "X：模拟时间步（整数）。", "Start_global": "Y：本步首次残差 / 固定 R_ref，无量纲。",
                   "Global_tolerance": f"全局非线性容差，恒为 {tol:g}。", "Time_s": "对应模拟时间，单位 s。"}
    book.table("Plot_Data", list(definitions), [[r["step"], r["start_global"], tol, r["time_s"]] for r in steps], definitions)
    chart = book.chart("Simulation time step", "Step-start residual / fixed R_ref", ymin=1e-15, ymax=20, xmax=72)
    book.series(chart, 1, "First logged residual in each step", ORANGE,
                marker={"type": "circle", "size": 3, "fill": {"color": ORANGE}, "border": {"color": ORANGE}})
    book.series(chart, 2, "Global nonlinear tolerance", GRAY, dashed=True)
    book.notes([("画图选列", "A 列作 X，B 列作 Y；C 列为水平容差线。使用 XY 散点带直线；Y 轴设为以 10 为底的对数轴。"),
                ("建议轴范围", "X：0–72；Y：1e-15–20。无需先对 B 列取 log10。"), *nonlinear_notes], definitions)
    results.append(book.close_and_verify())

    # B: already log-transformed reduction; its axis must stay linear.
    book = Book("B_Within_Step_Reduction.xlsx", "(B) Nonlinear reduction within each step",
                "residual_convergence (B) / residual_nonlinear 右图", source_hashes)
    definitions = {"Step": "X：模拟时间步。", "Orders_reduced": "Y：log10(Start_norm / Last_norm)，数量级降幅，已取过对数。",
                   "Nonlinear_iterations": "本步实际非线性迭代次数。", "Last_over_start": "本步最后残差 / 本步最初残差，无量纲。",
                   "Start_global": "本步首次残差 / 固定 R_ref。", "Last_global": "本步最后记录残差 / 固定 R_ref。",
                   "Iteration_change": "迭代次数变化的说明；仅第 4、23 步。"}
    rows = [[r["step"], r["reduction_orders"], r["iterations"], r["end_over_start"], r["start_global"], r["end_global"],
             {4: "4 to 3 iterations", 23: "3 to 2 iterations"}.get(r["step"])] for r in steps]
    book.table("Plot_Data", list(definitions), rows, definitions, decimal_columns=(1,))
    chart = book.chart("Simulation time step", "Orders reduced: log10(Start / Last)", kind="column", logarithmic=False, ymin=0, ymax=15.4)
    chart.add_series({"name": "Within-step reduction", "categories": ["Plot_Data", 1, 0, 71, 0],
                      "values": ["Plot_Data", 1, 1, 71, 1], "fill": {"color": BLUE}, "border": {"none": True},
                      "points": [{"fill": {"color": ORANGE}} if r["step"] in (4, 23) else {} for r in steps], "gap": 28})
    chart.set_legend({"none": True})
    book.raw("Source_Data", steps)
    book.notes([("画图选列", "A 列作类别 X，B 列作 Y，插入柱状图。Y 轴使用普通线性轴，不要再次设为对数轴。"),
                ("容差说明", "降幅为 10 表示本步残差缩小 10^10 倍；它不是统一的停止阈值，不应添加 Y=10 的停止线。"), *nonlinear_notes], definitions)
    results.append(book.close_and_verify())

    # C: two marker groups, with genuinely missing cells for the other group.
    book = Book("C_Linear_Solve_Termination.xlsx", "(C) Linear solve termination",
                "residual_convergence (C) / residual_linear_termination", source_hashes)
    definitions = {"Linear_solve_index": "X：全程线性求解序号。", "RTOL_true_relative": "Y1：相对容差停止组的最后真实相对残差。",
                   "ATOL_true_relative": "Y2：绝对容差停止组的最后真实相对残差。", "Relative_tolerance": "线性相对容差（水平线）。",
                   "All_true_relative": "全部求解的最后真实相对残差，便于检查。", "Stop_reason": "PETSc 实际终止原因。",
                   "Step": "所属模拟时间步。", "Nonlinear_iteration": "所属非线性迭代次数。"}
    rows = [[r["linear_solve_index"], r["final_true_relative_residual"] if r["reason"] == "CONVERGED_RTOL" else None,
             r["final_true_relative_residual"] if r["reason"] == "CONVERGED_ATOL" else None, rtol,
             r["final_true_relative_residual"], r["reason"], r["step"], r["nonlinear_iteration"]] for r in solves]
    book.table("Plot_Data", list(definitions), rows, definitions)
    chart = book.chart("Linear-solve index", "Final true relative residual", ymin=7e-11, ymax=5e-8, xmax=168)
    for col, name, color, symbol in [(1, "Relative-tolerance stop", "#176b91", "circle"),
                                     (2, "Absolute-tolerance stop", "#c55c26", "triangle")]:
        book.series(chart, col, name, color, connect=False,
                    marker={"type": symbol, "size": 4, "fill": {"color": color}, "border": {"color": color}})
    book.series(chart, 3, "Linear relative tolerance", GRAY, dashed=True)
    book.raw("Source_Data", solves)
    book.notes([("画图选列", "A 列为共同 X，B/C 列分别为两组 Y，仅画散点，不连接各次求解；D 列作容差线。Y 轴设为 log10。"),
                ("真实相对残差", "||b - A x|| / ||b||；取每次线性求解最后一次 GMRES 监视记录。"),
                ("线性停止条件", f"相对容差 {rtol:g} 或绝对容差 {atol:g}。绝对容差停止的点可以高于相对容差线，不能据此判为未收敛。"),
                ("绝对容差", "绝对容差针对内部缩放系统的绝对残差范数，不应直接画在相对残差的 Y 轴上。")], definitions)
    results.append(book.close_and_verify())

    # D: align by local GMRES iteration, leaving shorter-history tails empty.
    labels = [f"Step {solves[k - 1]['step']}, iteration {solves[k - 1]['nonlinear_iteration']}" for k in selected]
    book = Book("D_Selected_GMRES_Histories.xlsx", "(D) Selected GMRES histories",
                "residual_convergence (D) / residual_gmres_histories", source_hashes)
    definitions = {"GMRES_iteration": "X：各自线性求解内部的 GMRES 迭代编号，从 0 开始。"}
    for index, label in zip(selected, labels):
        definitions[f"Solve_{index}_true_relative"] = f"Y：{label}（线性求解序号 {index}）的真实相对残差。"
    definitions["Relative_tolerance"] = "线性相对容差，水平线。"
    size = max(len(v) for v in histories.values())
    rows = [[i, *[histories[k][i]["true_relative_residual"] if i < len(histories[k]) else None for k in selected], rtol] for i in range(size)]
    book.table("Plot_Data", list(definitions), rows, definitions)
    chart = book.chart("GMRES iteration", "True relative residual", xmax=720)
    for column, (index, label, color) in enumerate(zip(selected, labels, [BLUE, "#379778", "#c55c26"]), 1):
        book.series(chart, column, label, color, last=len(histories[index]))
    book.series(chart, 4, "Linear relative tolerance", GRAY, dashed=True)
    for index in selected:
        book.raw(f"Raw_Solve_{index}", histories[index])
    book.notes([("画图选列", "A 列为共同 X，B/C/D 列是三条独立曲线；E 列为容差线。插入 XY 散点带直线，Y 轴设为 log10。"),
                ("曲线长度", "；".join(f"{label}: {len(histories[index])} 点，迭代 0–{len(histories[index])-1}" for index, label in zip(selected, labels))),
                ("空白处理", "每条曲线结束后保持空白。不要补零，不要将三条曲线首尾相连，也不要替换为预条件残差。"),
                ("终止原因", f"第 167 次线性求解达到绝对容差 {atol:g} 后终止，因此曲线终点可高于相对容差线 {rtol:g}。")], definitions)
    results.append(book.close_and_verify())

    # E: a separate XY segment per step, never an inter-step last-value envelope.
    book = Book("E_Paired_Start_Last_Residuals.xlsx", "First and last logged residuals within each step",
                "residual_nonlinear_pairs 补充图", source_hashes)
    definitions = {"Step": "X：模拟时间步。", "Start_global": "Y1：首次记录残差 / 固定 R_ref。",
                   "Last_global": "Y2：最后记录残差 / 固定 R_ref。", "Global_tolerance": "全局非线性容差。",
                   "Nonlinear_iterations": "本步非线性迭代次数。"}
    book.table("Plot_Data", list(definitions), [[r["step"], r["start_global"], r["end_global"], tol, r["iterations"]] for r in steps], definitions)
    connectors = []
    for r in steps:
        connectors.extend([[r["step"], r["start_global"]], [r["step"], r["end_global"]], [None, None]])
    book.table("Connector_Data", ["Step", "Global_residual"], connectors)
    chart = book.chart("Simulation time step", "Logged residual / fixed R_ref", ymin=1e-15, ymax=20, xmax=72)
    for i, r in enumerate(steps):
        book.series(chart, 1, f"Step {r['step']} pair", "#b6bdc4", sheet="Connector_Data", first=3*i+1, last=3*i+2)
    book.series(chart, 1, "First logged record", ORANGE, connect=False,
                marker={"type": "circle", "size": 4, "fill": {"color": "white"}, "border": {"color": ORANGE}})
    book.series(chart, 2, "Last logged record", BLUE, connect=False,
                marker={"type": "circle", "size": 4, "fill": {"color": BLUE}, "border": {"color": BLUE}})
    book.series(chart, 3, "Global nonlinear tolerance", GRAY, dashed=True)
    chart.set_legend({"position": "bottom", "delete_series": list(range(len(steps))), "font": {"name": "Arial", "size": 11}})
    book.notes([("画图选列", "A 列为 X，B/C 列分别为首末残差散点，D 列为容差线；Y 轴设为 log10。"),
                ("配对竖线", "Connector_Data 每 2 点为同一步的竖线，中间以空行隔开。Excel 示例使用 71 个独立线段，避免跨时间步连接。"),
                ("禁止误连", "不要连接不同时间步的末次残差形成包络线；每条灰色竖线只表示同一步的首末残差变化。"), *nonlinear_notes], definitions)
    results.append(book.close_and_verify())

    assert {name: sha(OUT / name) for name in source_names} == source_hashes
    readme = """# 残差子图 Excel 数据

每幅图独立一个 .xlsx；第一张表 Plot_Data 可直接选列绘图，Example_Chart 是可编辑的 Excel 示例图，Notes 提供中文定义与操作说明。
所有数据来自现有审核后的 CSV，无需重新运行流场。示例图复用同一批数据，未复制原 PNG 中全部文字注释。

| 文件 | 对应子图 | X | Y | Y 轴 |
|---|---|---|---|---|
| A_Global_Transient_Evolution.xlsx | A 全局瞬态演化 | A 列 Step | B 列 Start_global；C 列容差 | log10 |
| B_Within_Step_Reduction.xlsx | B 每步非线性降幅 | A 列 Step | B 列 Orders_reduced | 线性（数据已取 log10） |
| C_Linear_Solve_Termination.xlsx | C 线性求解终止 | A 列 Linear_solve_index | B/C 列分别为相对/绝对容差停止组；D 列容差 | log10 |
| D_Selected_GMRES_Histories.xlsx | D 代表性 GMRES 历史 | A 列 GMRES_iteration | B/C/D 列三条独立曲线；E 列容差 | log10 |
| E_Paired_Start_Last_Residuals.xlsx | 首末残差配对补充图 | A 列 Step | B/C 列首末残差；D 列容差 | log10 |

A、B、E 各 71 个时间步；C 共 167 次线性求解；D 三条曲线分别有 704、673、279 个点。
所有残差列保存为数值，可以直接用于 Excel、Origin 或 Python；不能把分组空值或曲线末端空值改成 0。
对数轴直接使用残差原值，不需预先取 log10；B 图例外，其列值已是数量级降幅，使用线性轴。
配对图每条竖线只连接同一步的两个点，不跨步连线。归一化参考值全程固定。
绝对容差终止允许最终相对残差高于相对容差线。残差大小不等于物理解误差。

验证：EXCEL_VALIDATION.json 记录输入/输出 SHA256、数据行数、数值单元格回读和轴类型检查。
原始 CSV 保留原始有效位，Excel 数值约有 15 位有效数字；无需将科学计数法数字转换成文本。
"""
    (DEST / "README_ZH.md").write_text(readme, encoding="utf-8")
    validation = {"all_pass": True, "created_utc": datetime.now(timezone.utc).isoformat(),
                  "source_directory": str(OUT), "source_sha256": source_hashes,
                  "source_unchanged": True, "exporter_sha256": sha(Path(__file__)),
                  "xlsxwriter_version": xlsxwriter.__version__, "workbooks": results,
                  "linear_stop_counts": dict(Counter(r["reason"] for r in solves)),
                  "selected_histories_points": {str(k): len(v) for k, v in histories.items()},
                  "gmres_final_points_match_linear_summary": True,
                  "excel_charts": "Native editable charts; data, ZIP/XML, series count and log axes checked, not desktop-rendered."}
    (DEST / "EXCEL_VALIDATION.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    files = [DEST / r["file"] for r in results] + [DEST / "README_ZH.md", DEST / "EXCEL_VALIDATION.json"]
    (DEST / "SHA256SUMS.txt").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    archive_path = DEST / "Residual_Plot_Data_Excel.zip"
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        for path in files + [DEST / "SHA256SUMS.txt"]:
            archive.write(path, arcname=path.name)
    with ZipFile(archive_path) as archive:
        assert archive.testzip() is None
        assert len(archive.namelist()) == 8
    print(json.dumps({"archive": str(archive_path), "all_pass": True, "workbooks": results}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
