"""Sequential follow-up jobs keep a single resource owner."""
from .parallel import main as parallel
from .extension import main as extension
from .profile_extra import main as profile
if __name__=='__main__':
    parallel()
    extension()
    profile()
