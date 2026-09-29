"""Compatibility module routing every new production run to the fixed-window CLI.

Legacy --minimum-exits, --post-infusion and --fixed-infusion are rejected.
Historical coverage runs remain reproducible with their recorded code commit.
"""
from .run_fixed_window_population import main
if __name__=='__main__':main()
