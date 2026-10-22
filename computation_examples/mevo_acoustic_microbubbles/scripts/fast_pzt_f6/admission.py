"""Run final software admission sequentially on the sole GPU."""

import argparse
import sys
from . import checks, efficiency


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    label = parser.parse_args().label
    for module, suffix in [(checks, "_checks"), (efficiency, "_efficiency")]:
        sys.argv = [module.__name__, "--label", label + suffix]
        module.main()


if __name__ == "__main__":
    main()
