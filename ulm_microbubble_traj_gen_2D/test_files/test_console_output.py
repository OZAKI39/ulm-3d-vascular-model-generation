from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import unittest
from unittest.mock import patch

from ulm_microbubble_traj_gen_2D.utils.runtime.console_output import (
    print_banner,
    print_command,
    print_error,
    print_key_values,
    print_stage,
    print_warning,
)


class ConsoleOutputTests(unittest.TestCase):
    def test_banner_has_stable_rules_title_and_subtitle(self) -> None:
        stream = StringIO()
        with patch(
            "ulm_microbubble_traj_gen.utils.runtime.console_output._console_width",
            return_value=72,
        ):
            print_banner("Trajectory generation", subtitle="Mode: quick test", file=stream)

        lines = stream.getvalue().splitlines()
        self.assertEqual(lines[0], "=" * 72)
        self.assertEqual(lines[1], " Trajectory generation")
        self.assertEqual(lines[2], " Mode: quick test")
        self.assertEqual(lines[3], "=" * 72)

    def test_stage_uses_zero_padded_progress_index(self) -> None:
        stream = StringIO()
        print_stage(3, 8, "Build domain", detail="128 x 96", file=stream)

        output = stream.getvalue()
        self.assertIn("[03/08] Build domain", output)
        self.assertIn("128 x 96", output)

    def test_key_values_align_colons_and_indent_multiline_values(self) -> None:
        stream = StringIO()
        print_key_values(
            [("Short", 1), ("Longer label", "first\nsecond")],
            file=stream,
        )

        lines = stream.getvalue().splitlines()
        self.assertEqual(lines[0].index(":"), lines[1].index(":"))
        self.assertTrue(lines[2].endswith("second"))
        self.assertGreater(lines[2].index("second"), lines[1].index(":"))

    def test_command_is_copy_ready_without_a_key_value_separator(self) -> None:
        stream = StringIO()
        print_command(r"D:\anaconda3\envs\pmp\python.exe -m package.viewer", file=stream)

        self.assertEqual(
            stream.getvalue(),
            "  D:\\anaconda3\\envs\\pmp\\python.exe -m package.viewer\n",
        )

    def test_warning_wraps_with_one_visible_prefix(self) -> None:
        stream = StringIO()
        with patch(
            "ulm_microbubble_traj_gen.utils.runtime.console_output._console_width",
            return_value=72,
        ):
            print_warning("word " * 40, file=stream)

        output = stream.getvalue()
        self.assertEqual(output.count("[WARNING]"), 1)
        self.assertGreater(len(output.splitlines()), 1)

    def test_error_defaults_to_stderr(self) -> None:
        stdout = StringIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            print_error("generation stopped")

        self.assertEqual(stdout.getvalue(), "")
        self.assertIn("[ERROR] generation stopped", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
