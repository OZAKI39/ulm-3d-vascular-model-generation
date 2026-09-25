# 浓度语义追溯：total / conditional 已解决，稳态绝对值为模型假设

旧 `PARTICLE7_INLET_POPULATION_V0.json` 将 8.5e12 m⁻³ 明确命名为 `NOMINAL_WELL_MIXED_POST_BOLUS_ANCHOR`。原论文的动物实验部分给出 32.6 g 平均体重、2×10⁷ 个 microbubbles 的 bolus，未将该剂量限定为 D≤4 µm。该文没有测量本血管的连续灌注稳态浓度。[原论文 2.6.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12445601/)

血容量取历史合同的 72 mL/kg；独立机构指南也给出小鼠这一估计值。它是估算而非本实验逐鼠测量。[University of Michigan guideline](https://az.research.umich.edu/animalcare/guidelines/guidelines-blood-collection)

计算：V_blood=0.0326×72=2.3472 mL，C_total_nominal=2e7/2.3472×1e6=8.52079072938e+12 m⁻³。继承原舍入锚点 C_total=8.5e12 m⁻³。不能把这个全部 SonoVue 浓度直接写成 conditional modeled 浓度。

冻结分布的 F_original(4 µm)=0.989750470288813。所以 **C_modeled_D≤4=8412878997454.91 m⁻³**。源分布条件化与源强度扣除被排除尺寸质量必须同时发生。数学等价于从全体源删除 D>4 的 marks，再对真实入口做有限尺寸 thinning。

该文 bolus dose 支持“全部泡的名义剂量锚点”这一语义；本文进一步假设连续灌注把此名义浓度维持为常量。这一步是显式 `MODEL_ASSUMPTION`，不是从论文推得的生理 ground truth。未新增 pharmacokinetics、clearance、泡寿命或循环再入模型。若实际浓度变化，当前绝对时间需重新建模；固定稳态 C 的缩放不改变 normalized source/entering joint law，但按比例改变 lambda 和时间尺度。

100k 模型时钟约 8.85 天只用于统计精度，不表示 bolus 在活体持续这一时间。100k 入口生成无轨迹推进；30 smoke 为独立事件的局部 residence integration，非连续多日 many-body 模拟。

原始分布 SHA `2c9f783c6169421b06c057e16653878e7385139a179682c0648519e72a9f5198`；sampler SHA `c3f5ab75e20503a495cd107cef9b63240bda890e90f601d89b1f28fa8f292ece`；源合同 SHA `82ca08f51b3efef290c53beda53707b969b0e45e49ee6bb3bc636f3fb2f792d9`。原文全文通过 Europe PMC fullTextXML 核验，保存核验字段和 XML SHA，未将论文全文拷入仓库：[核验记录](data/concentration_literature_verification.json)。
