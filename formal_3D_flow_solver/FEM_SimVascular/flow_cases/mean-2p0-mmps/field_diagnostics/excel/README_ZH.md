# 残差子图 Excel 数据

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
