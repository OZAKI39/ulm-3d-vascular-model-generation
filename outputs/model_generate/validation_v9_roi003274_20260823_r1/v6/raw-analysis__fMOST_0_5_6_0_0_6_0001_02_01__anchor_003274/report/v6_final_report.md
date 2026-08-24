# v6 Continuous-Field Transition 最终核验

总体状态：**PASS**

## 方法约束

- source SWC、source radius、CORE/context 拓扑和 CUT_PORT 映射保持不变。
- Junction 使用单个 continuous implicit field 与一次 marching cubes。
- surface loop stitch 仅作为 v5 regression，不参与 v6 正式表面。
- merge 仅位于 w=1 的 PURE_BRANCH overlap；未做全局 smoothing。

## 关键结果

- Junction interface P99: 13.2952° → 9.1226°。
- Junction interface max: 17.2004° → 12.9996°。
- Port interface max: 7.69069° → 2.19069°。
- Silhouette curvature variation mean: 0.259068 → 0.259017。
- Topology defects: boundary=0, nonmanifold=0, self-intersection=0, internal faces=0, internal caps=0, degenerate=0。

最终 PASS 不由测试用例单独决定；以上定量结果与 figures/ 中同相机图共同构成判据。
