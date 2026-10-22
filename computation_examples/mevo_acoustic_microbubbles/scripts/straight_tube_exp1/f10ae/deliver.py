"""Final evidence assembly after frozen scientific attempts; no new simulation."""
from .common import BASE, read

def main():
    assert (BASE/'LOW_DOMAIN_FOLLOWUP.json').exists()
    from . import interpret, figures, finalize, report
    interpret.main()
    figures.main()
    finalize.main()
    report.main()
    finalize.main()
    print('DELIVERY_ASSEMBLED', read(BASE/'INTERPRETATION.json')['physical_attempts'], flush=True)

if __name__ == '__main__':
    main()
