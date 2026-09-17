"""
Tinit 子树选择、独立阶段树复制，以及阶段结果合并。

分阶段生长的难点在于：某个阶段可能只在全局树的一部分上继续长。
为了不直接破坏全局树，流程是：
1. 从全局树选择 Tinit 子集；
2. 把这个子集复制成一棵独立的 `DCCOTree`；
3. 在这棵阶段树上运行普通生长算法；
4. 把阶段树中新长出的部分合并回全局树。

`stage_to_global` / `global_to_stage` 这两个字典就是为了在阶段树 ID 和全局 ID
之间翻译，避免合并时接错父子关系。
"""

from dataclasses import dataclass, field

import numpy as np

from ..core.models import Vessel
from ..core.tree import DCCOTree


@dataclass
class StageRegistry:
    """记录每个命名 stage 在全局树中修改或新增了哪些 vessel。"""

    # key 是 stage 名，value 是该 stage 结束后产生或修改的全局 vessel id 集合。
    stage_outputs: dict = field(default_factory=dict)
    # previous_stage 记录最近完成的阶段名，支持 tinit: previous 写法。
    previous_stage: str = None


@dataclass(frozen=True)
class StageInitialState:
    """Initial tree and ID mappings prepared for one growth stage."""
    """一个 stage 的独立阶段树，以及阶段树 ID 到全局 ID 的对应关系。"""

    # 独立 DCCOTree，后续普通生长算法只作用在这棵树上。
    tree: object
    # 阶段树中的 vessel id -> 全局树中的 vessel id。
    stage_to_global: dict
    # 用户真正选择为可继续分支的全局 vessel id。
    selected_global_ids: set
    # 阶段树初始包含的所有全局 vessel id，包含必要祖先/后代。
    initial_global_ids: set


@dataclass(frozen=True)
class StageMergeResult:
    """阶段树合并回全局树后的结果统计。"""

    # 当前 stage 对外暴露的输出 id，供后续 stage 继续引用。
    stage_output_ids: set
    # 阶段树中新长出的、追加到全局树的 vessel id。
    new_global_ids: set
    # 原本已存在但被当前 stage 改写几何/血流/连接关系的全局 vessel id。
    modified_global_ids: set


def _topological_subset_order(vessels, ids):
    """返回某个子集内部的拓扑顺序，确保复制成新树时父节点先出现。"""
    # 空子集没有顺序，直接返回空列表。
    if not ids:
        return []
    # 子集内部的根是：父节点不存在，或父节点不在子集里。
    roots = sorted(
        vid for vid in ids
        if vessels[vid].parent_id < 0 or vessels[vid].parent_id not in ids
    )
    # order 保存最终拓扑序，seen 防止重复访问。
    order = []
    seen = set()
    # stack 用 reversed 初始化，是为了 pop 时仍按升序访问 root。
    stack = list(reversed(roots))
    # 深度优先遍历子集。
    while stack:
        vid = stack.pop()
        if vid in seen:
            continue
        # 记录访问状态和输出顺序。
        seen.add(vid)
        order.append(vid)
        # 只遍历同样在子集里的子节点。
        child_ids = sorted(c for c in vessels[vid].children if c in ids)
        # reversed 后压栈，保证弹出时仍是升序。
        stack.extend(reversed(child_ids))
    # 返回父先于子的顺序，复制成新树时才能正确设置 parent_id。
    return order


def select_tinit_vessels(vessels, spec, registry):
    """
    根据 TinitSpec 选择初始子树和可继续分叉的 vessel。

    `initial_ids` 是 stage 阶段树里会出现的所有 vessel；
    `branchable_ids` 是其中允许继续长新分支的 vessel。包含祖先但不允许
    祖先分叉，可以保持树连通，同时避免无关主干继续生长。
    """
    # empty 模式表示当前阶段从空树开始，不继承任何全局 vessel。
    if spec.mode == "empty":
        selected = set()
    # global 模式选择当前全局树中的所有 vessel。
    elif spec.mode == "global":
        selected = {v.vid for v in vessels}
    elif spec.mode in {"previous", "stages"}:
        if spec.mode == "previous":
            names = (registry.previous_stage,)
        else:
            names = spec.stage_names
        selected = set()
        for name in names:
            if name not in registry.stage_outputs:
                raise ValueError(f"Unknown Tinit stage: {name!r}")
            selected.update(registry.stage_outputs[name])
    elif spec.mode == "ids":
        selected = set(spec.vessel_ids)
    elif spec.mode == "terminals":
        selected = set(spec.terminal_ids)
    elif spec.mode == "query":
        selected = {v.vid for v in vessels}
    else:
        raise ValueError(f"Unsupported Tinit mode: {spec.mode!r}")

    if spec.mode not in {"ids", "terminals"}:
        selected.update(spec.vessel_ids)
        terminal_ids = set(spec.terminal_ids)
        selected.update(terminal_ids)
    elif spec.mode == "terminals":
        terminal_ids = selected
    else:
        terminal_ids = set()

    unknown_ids = sorted(vid for vid in selected if vid < 0 or vid >= len(vessels))
    if unknown_ids:
        raise ValueError(f"Tinit references unknown vessel IDs: {unknown_ids}")

    non_terminals = sorted(vid for vid in terminal_ids if vessels[vid].children)
    if non_terminals:
        raise ValueError(
            f"Tinit terminal IDs are not terminal vessels: {non_terminals}"
        )

    # roles 过滤允许只保留某些功能角色的 vessel。
    if spec.roles:
        roles = set(spec.roles)
        selected = {vid for vid in selected if vessels[vid].role in roles}
    # branching_modes 过滤允许只保留某些生成来源的 vessel。
    if spec.branching_modes:
        modes = set(spec.branching_modes)
        selected = {vid for vid in selected if vessels[vid].branching_mode in modes}

    # initial_ids 是阶段树中实际需要包含的 vessel 集合。
    initial_ids = set(selected)
    # 包含祖先可保证选中子树和原全局根之间不断开。
    if spec.include_ancestors:
        for start in selected:
            vid = vessels[start].parent_id
            while vid >= 0 and vid not in initial_ids:
                initial_ids.add(vid)
                vid = vessels[vid].parent_id

    # 包含后代可完整保留选中局部子树的下游结构。
    if spec.include_descendants:
        stack = list(selected)
        seen = set()
        while stack:
            vid = stack.pop()
            if vid in seen:
                continue
            seen.add(vid)
            initial_ids.add(vid)
            stack.extend(vessels[vid].children)

    # global 或 branch_from_all 模式下，阶段树中的每段都可作为新分叉候选。
    if spec.mode == "global" or spec.branch_from_all:
        branchable_ids = set(initial_ids)
    # 其他模式只允许原始 selected 集合继续分叉，祖先/后代只是为了连通和完整。
    else:
        branchable_ids = set(selected) & initial_ids

    # 返回阶段树需要包含的集合，以及允许继续生长的集合。
    return initial_ids, branchable_ids


def _copy_vessel(v, vid, parent_id, children):
    """复制一个 Vessel，同时替换它在目标树中的 id、父 id 和 children。"""
    # copy() 用于复制坐标数组，避免阶段树和全局树共享同一块可变内存。
    return Vessel(
        vid=vid,
        parent_id=parent_id,
        children=children,
        x_p=v.x_p.copy(),
        x_d=v.x_d.copy(),
        radius=v.radius,
        branching_mode=v.branching_mode,
        role=v.role,
        is_main_trunk=v.is_main_trunk,
        prescribed_outflow=v.prescribed_outflow,
        flow_rate=v.flow_rate,
        mean_velocity=v.mean_velocity,
        flow_conservation_residual=v.flow_conservation_residual,
        murray_residual=v.murray_residual,
    )


def build_stage_initial_state(cfg, global_vessels, tinit, registry):
    """Build the isolated initial tree and ID mappings for one growth stage."""
    """把 stage 指定的 Tinit 子集复制成独立 DCCOTree。"""
    # 根据 tinit 规则得到“要复制哪些 vessel”和“哪些 vessel 允许继续分叉”。
    initial_ids, branchable_ids = select_tinit_vessels(global_vessels, tinit, registry)
    # 创建一棵空的阶段树，配置使用当前阶段自己的 cfg。
    tree = DCCOTree(cfg)
    # 如果初始集合为空，直接返回空树，后续生长会按普通空树逻辑初始化。
    if not initial_ids:
        return StageInitialState(
            tree=tree,
            stage_to_global={},
            selected_global_ids=set(),
            initial_global_ids=set(),
        )

    # 拓扑顺序保证复制父节点时，它的父 vessel 已经有阶段树 id。
    order = _topological_subset_order(global_vessels, initial_ids)
    # global_to_stage 方便从全局 id 找阶段树 id。
    global_to_stage = {}
    # stage_to_global 方便合并时从阶段树 id 找回全局 id。
    stage_to_global = {}
    # 第一遍：复制 vessel 本体，并建立 id 映射。
    for gid in order:
        # src 是全局树中的原始 vessel。
        src = global_vessels[gid]
        # 阶段树 id 按追加顺序生成。
        stage_id = len(tree.vessels)
        # 双向映射必须同时记录，后续重连 parent/children 都会用到。
        global_to_stage[gid] = stage_id
        stage_to_global[stage_id] = gid
        # 如果父节点不在复制出的子集内，则父 id 设为 -1，形成该子树局部根。
        if src.parent_id in global_to_stage:
            parent_stage = global_to_stage[src.parent_id]
        else:
            parent_stage = -1
        # children 暂时传空列表，第二遍再统一填充。
        stage_vessel = _copy_vessel(src, stage_id, parent_stage, [])
        # 不允许继续分叉的祖先/后代保留在树里，但标记为 non_branching。
        if gid not in branchable_ids:
            stage_vessel.branching_mode = "non_branching"
        # 把复制出的 vessel 加入阶段树。
        tree.vessels.append(stage_vessel)

    # 第二遍：根据全局 children 和 id 对应关系重建阶段树的 children。
    for gid in order:
        # 找到全局 id 对应的阶段树 id。
        stage_id = global_to_stage[gid]
        # 只保留同样被复制进阶段树的子节点。
        tree.vessels[stage_id].children = [
            global_to_stage[child_id]
            for child_id in global_vessels[gid].children
            if child_id in global_to_stage
        ]

    # 阶段树准备完成后清空缓存，避免旧拓扑查询结果被复用。
    tree._invalidate_tree_cache()
    # 根据当前 cfg 的半径模式重新同步半径，确保阶段树内部状态一致。
    tree._apply_radii()
    # 返回阶段树和映射信息，供阶段生长与合并使用。
    return StageInitialState(
        tree=tree,
        stage_to_global=stage_to_global,
        selected_global_ids=set(branchable_ids),
        initial_global_ids=set(initial_ids),
    )


def _copy_stage_geometry(dst, src, keep_branching_metadata):
    """把 stage 阶段树中的几何和血流结果复制回全局 vessel。"""
    # 几何端点必须复制数组，避免全局树和阶段树共享同一数组。
    dst.x_p = src.x_p.copy()
    dst.x_d = src.x_d.copy()
    # 半径和血流字段代表阶段优化后的结果，直接覆盖全局值。
    dst.radius = src.radius
    dst.prescribed_outflow = src.prescribed_outflow
    dst.flow_rate = src.flow_rate
    dst.mean_velocity = src.mean_velocity
    dst.flow_conservation_residual = src.flow_conservation_residual
    dst.murray_residual = src.murray_residual
    dst.is_main_trunk = src.is_main_trunk
    # 只有新 vessel 才保留 stage 内部元数据；旧 vessel 的角色通常来自全局树。
    if keep_branching_metadata:
        dst.branching_mode = src.branching_mode
        dst.role = src.role


def merge_stage_tree(global_vessels, stage_vessels, stage_to_global):
    """
    把完成生长的 stage 阶段树合并回全局树。

    这个函数要同时处理三类情况：
    - 阶段树中的新 vessel：追加到全局树；
    - 阶段树中对应全局旧 vessel 的节点：更新几何/半径/血流；
    - 旧 vessel 被 stage 分裂：把原本挂在旧远端的全局孩子重新接到
      successor 段，保持全局树拓扑正确。
    """
    # 记录合并前旧全局段的远端；如果旧段被分裂，后面要靠它找到 successor。
    old_distal = {gid: global_vessels[gid].x_d.copy() for gid in stage_to_global.values()}
    # 记录合并前旧全局段的孩子；合并后可能需要把这些孩子重新挂接。
    old_children = {gid: list(global_vessels[gid].children) for gid in stage_to_global.values()}
    # full_map 最终包含所有阶段树 id 到全局树 id 的对应关系，包括新 vessel。
    full_map = dict(stage_to_global)
    # new_global_ids 记录本次新增到全局树的 vessel。
    new_global_ids = set()
    # modified_global_ids 记录本次修改过的已有全局 vessel。
    modified_global_ids = set()

    # 第一轮：先把阶段树里还没有全局 id 的新 vessel 追加进全局树。
    stage_order = _topological_subset_order(stage_vessels, {v.vid for v in stage_vessels})
    for stage_id in stage_order:
        # 已在 stage_to_global 中的 vessel 是旧全局段，不需要新建。
        if stage_id in full_map:
            continue
        # stage_vessel 是阶段树中新长出来的 vessel。
        stage_vessel = stage_vessels[stage_id]
        # 新全局 id 采用追加到列表末尾的方式生成。
        new_id = len(global_vessels)
        # 记录新阶段树 id 到新全局树 id 的对应关系。
        full_map[stage_id] = new_id
        # 新 vessel 的父节点也要从阶段树 id 转成全局树 id。
        if stage_vessel.parent_id < 0:
            parent_id = -1
        else:
            parent_id = full_map[stage_vessel.parent_id]
        # 复制新 vessel 到全局树，children 稍后统一重建。
        global_vessels.append(_copy_vessel(stage_vessel, new_id, parent_id, []))
        # 记录新增 id，供日志和后续 stage 引用。
        new_global_ids.add(new_id)

    # 第二轮：把阶段树中所有 vessel 的几何、血流和 parent 关系写回全局树。
    for stage_vessel in stage_vessels:
        # 找到当前阶段树 vessel 对应的全局 vessel。
        gid = full_map[stage_vessel.vid]
        dst = global_vessels[gid]
        # 保存写回前的端点和 children，用于判断是否真的改动了旧 vessel。
        before_xp = dst.x_p.copy()
        before_xd = dst.x_d.copy()
        before_children = set(dst.children)
        # 新 vessel 可以继承阶段元数据；旧 vessel 通常只更新几何和血流。
        keep_branching_metadata = stage_vessel.vid not in stage_to_global
        # 写回几何、半径和血流字段。
        _copy_stage_geometry(dst, stage_vessel, keep_branching_metadata)
        # 阶段树中有父节点时，需要把父 id 翻译成全局 id。
        if stage_vessel.parent_id >= 0:
            dst.parent_id = full_map[stage_vessel.parent_id]
        # 新的根节点保持 parent_id=-1。
        elif stage_vessel.vid not in stage_to_global:
            dst.parent_id = -1

        # mapped_children 是阶段树 children 转成全局 id 后的集合。
        mapped_children = [full_map[child_id] for child_id in stage_vessel.children]
        # 新增加的阶段树 child 说明旧 vessel 在本 stage 中被接上了新分支。
        added_stage_children = set(mapped_children) - before_children
        # 如果旧 vessel 的端点或孩子关系发生变化，就把它记为 modified。
        if (not np.allclose(before_xp, dst.x_p, atol=1e-7, rtol=0.0)
                or not np.allclose(before_xd, dst.x_d, atol=1e-7, rtol=0.0)
                or added_stage_children):
            if stage_vessel.vid in stage_to_global:
                modified_global_ids.add(gid)

    # 第三轮：处理“旧全局段被 stage 分裂”时的旧孩子重接问题。
    for stage_id, gid in stage_to_global.items():
        # stage_vessel 是旧全局段在 stage 树中的对应段。
        stage_vessel = stage_vessels[stage_id]
        # 将阶段树当前孩子转成全局 id，便于和旧 children 比较。
        mapped_children = {full_map[child_id] for child_id in stage_vessel.children}
        # 如果旧远端没有移动，说明没有发生分裂，不需要重接。
        if np.allclose(
            old_distal[gid],
            global_vessels[gid].x_d,
            atol=1e-7,
            rtol=0.0,
        ):
            continue

        successor = None
        parent = global_vessels[gid]
        for child_id in mapped_children:
            child = global_vessels[child_id]
            if (
                np.allclose(child.x_p, parent.x_d, atol=1e-7, rtol=0.0)
                and np.allclose(
                    child.x_d,
                    old_distal[gid],
                    atol=1e-7,
                    rtol=0.0,
                )
            ):
                successor = child_id
                break

        # 找不到 successor 时保守跳过，避免错误改接拓扑。
        if successor is None:
            continue
        # 遍历合并前旧段的孩子，判断哪些需要从旧段改挂到 successor。
        for child_id in old_children[gid]:
            # 已经属于 stage 新拓扑的孩子不再处理。
            if child_id in mapped_children:
                continue
            # 只有原本接在旧远端的孩子才应该改接到 successor。
            if (
                global_vessels[child_id].parent_id == gid
                and np.allclose(
                    global_vessels[child_id].x_p,
                    old_distal[gid],
                    atol=1e-7,
                    rtol=0.0,
                )
            ):
                global_vessels[child_id].parent_id = successor
                modified_global_ids.add(child_id)

    for vessel in global_vessels:
        vessel.children = []
    for vessel in global_vessels:
        if vessel.parent_id >= 0:
            global_vessels[vessel.parent_id].children.append(vessel.vid)
    # stage 输出包括新增 vessel 和被修改的旧 vessel。
    stage_output_ids = set(new_global_ids) | set(modified_global_ids)
    # 返回合并结果统计。
    return StageMergeResult(
        stage_output_ids=stage_output_ids,
        new_global_ids=new_global_ids,
        modified_global_ids=modified_global_ids,
    )
