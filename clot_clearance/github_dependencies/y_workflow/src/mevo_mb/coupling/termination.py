"""Objective-driven termination; classification alone is not physical terminality."""
from dataclasses import dataclass


@dataclass(frozen=True)
class SimulationTerminationPolicy:
    objective: str = "finite_observation_with_contact_metrics"
    terminal_outlet_proven: bool = False
    no_return_proven: bool = False
    excluded_pair_influence_proven: bool = False
    require_downstream_metrics: bool = True

    @property
    def terminal_deactivation_safe(self):
        return (self.terminal_outlet_proven and self.no_return_proven
                and self.excluded_pair_influence_proven and not self.require_downstream_metrics)

    def reason(self, *, time_s, observation_limit_s, total_particles, terminally_assigned):
        if time_s >= observation_limit_s-1e-15:
            return "observation_limit_reached"
        if (self.objective == "terminal_branch_assignment_only" and self.terminal_deactivation_safe
                and total_particles > 0 and terminally_assigned == total_particles):
            return "all_particles_physically_terminal_and_objective_complete"
        return None
