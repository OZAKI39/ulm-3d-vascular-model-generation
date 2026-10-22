"""Final readable scientific labels, unchanged physical input trajectories."""
from . import videos, figures, report, finalize

def main():
    videos.main()
    figures.main()
    finalize.main()
    report.main()
    finalize.main()

if __name__ == '__main__':
    main()
