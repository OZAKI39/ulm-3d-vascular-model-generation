"""Default: configure PZT1 but leave it OFF. Pass --enable to emit for 2 s.

Use one control process at a time. Stop the GUI before running this example.
"""
import argparse
import time
from afg_control import Rig, ChannelSettings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--enable", action="store_true", help="Explicitly enable PZT1 at 0.1 Vpp for 2 seconds")
    args = parser.parse_args()
    with Rig() as rig:
        rig.apply({"PZT1": ChannelSettings(frequency_hz=2.5e6, amplitude_vpp=0.1, phase_deg=0)})
        if args.enable:
            rig.set_output("PZT1", True)
            time.sleep(2)
        print(rig.read_channel("PZT1"))
    # Context exit attempts to turn all four outputs OFF, including on exceptions.


if __name__ == "__main__":
    main()
