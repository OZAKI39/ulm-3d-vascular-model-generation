"""
Run the staged vascular-growth pipeline.
"""

from ..core.tree import DCCOTree
from .stage_union import (
    StageRegistry,
    build_stage_initial_state,
    merge_stage_tree,
)


def generate_staged_tree(stages):
    """
    Grow each configured stage and return the final global tree.
    """
    registry = StageRegistry()
    global_vessels      = []
    for stage in stages:
        # Copy the selected Tinit subtree into an isolated stage-growth state.
        initial_state = build_stage_initial_state(stage.cfg, global_vessels, stage.tinit, registry)

        stage_tree              = initial_state.tree
        stage_tree._fr_runtime  = 1.0

        # formally grow the stage tree to the number of terminals specified in the stage configuration.
        # start from root if tinit is None, otherwise start from the terminals of the provided initial tree.
        stage_tree.grow_additional_terminals(stage.cfg.n_terminals)

        result = merge_stage_tree(
            global_vessels,
            stage_tree.vessels,
            initial_state.stage_to_global,
        )

        global_tree         = DCCOTree(stage.cfg)
        global_tree.vessels = global_vessels
        
        global_tree._invalidate_tree_cache()
        global_tree._apply_radii()
        registry.stage_outputs[stage.name] = set(result.stage_output_ids)
        registry.previous_stage = stage.name

    return global_tree
