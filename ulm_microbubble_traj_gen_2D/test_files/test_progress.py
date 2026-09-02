from __future__ import annotations

import unittest
from unittest import mock

from ulm_microbubble_traj_gen_2D.utils.runtime import progress


class RuntimeProgressTests(unittest.TestCase):
    def test_stage_progress_creates_and_completes_live_current_module_bar(
        self,
    ) -> None:
        stage_bar = mock.Mock()
        current_bar = mock.Mock()
        with mock.patch.object(
            progress,
            "create_progress_bar",
            side_effect=[stage_bar, current_bar],
        ) as create:
            stage_progress = progress.create_stage_progress_bar(4, 10)
            stage_progress.set_postfix(module="Gmsh mesh", refresh=True)
            stage_progress.set_submodule_progress(
                completed=3,
                total=8,
                submodule="create boundary geometry",
            )
            stage_progress.update()
            stage_progress.close()

        self.assertEqual(create.call_count, 2)
        self.assertEqual(create.call_args_list[0].kwargs["position"], 0)
        self.assertEqual(create.call_args_list[1].kwargs["position"], 1)
        self.assertIn(
            "Gmsh mesh",
            create.call_args_list[1].kwargs["description"],
        )
        stage_bar.set_postfix.assert_called_once_with(
            module="Gmsh mesh",
            refresh=True,
        )
        self.assertEqual(current_bar.total, 8)
        self.assertEqual(
            current_bar.update.call_args_list,
            [mock.call(3), mock.call(5)],
        )
        current_bar.set_postfix.assert_called_once_with(
            submodule="create boundary geometry",
            refresh=True,
        )
        current_bar.close.assert_called_once_with()
        stage_bar.update.assert_called_once_with(1)
        stage_bar.close.assert_called_once_with()

    def test_starting_next_module_closes_incomplete_previous_bar(self) -> None:
        stage_bar = mock.Mock()
        first_bar = mock.Mock()
        second_bar = mock.Mock()
        with mock.patch.object(
            progress,
            "create_progress_bar",
            side_effect=[stage_bar, first_bar, second_bar],
        ):
            stage_progress = progress.create_stage_progress_bar(3, 2)
            stage_progress.set_postfix(module="continuous geometry")
            stage_progress.set_postfix(module="vessel rasterization")
            stage_progress.close()

        first_bar.update.assert_not_called()
        first_bar.close.assert_called_once_with()
        second_bar.update.assert_not_called()
        second_bar.close.assert_called_once_with()

    def test_particle_progress_is_the_primary_stage_six_bar(self) -> None:
        with mock.patch.object(
            progress,
            "create_progress_bar",
            return_value=mock.Mock(),
        ) as create:
            progress.create_particle_progress_bar(250)

        create.assert_called_once_with(
            250,
            description="Stage 06 simulation",
            unit="substep",
            position=0,
            leave=True,
        )


if __name__ == "__main__":
    unittest.main()
