"""Term-by-term SI dimension ledger (M,L,T exponents), not runtime units."""

from dataclasses import dataclass


@dataclass(frozen=True)
class D:
    m: int = 0
    length: int = 0
    t: int = 0

    def __mul__(self, o):
        return D(self.m + o.m, self.length + o.length, self.t + o.t)

    def __truediv__(self, o):
        return D(self.m - o.m, self.length - o.length, self.t - o.t)

    def __pow__(self, n):
        return D(self.m * n, self.length * n, self.t * n)


def test_zhang_every_radial_translational_shell_force_dimension():
    length, time, mass = D(length=1), D(t=1), D(m=1)
    rho = mass / length**3
    p = mass / length / time**2
    speed = length / time
    acc = length / time**2
    mu = mass / length / time
    sigma = chi = mass / time**2
    epsilon = length
    mus = mu
    radial = [
        length * acc,
        speed**2,
        p / rho,
        length / (length / time) * (p / time) / rho,
        length**2 * acc / length,
        length * speed**2 / length,
        length**2 * speed**2 / length**2,
        length**3 * acc / length**2,
        length**3 * speed**2 / length**3,
    ]
    assert all(v == D(length=2, t=-2) for v in radial)
    force = mass * length / time**2
    translation = [
        length * acc,
        speed**2,
        (length * length**2 * speed) / time / length**2,
        length**2 * length**2 * acc / length**3,
        length**3 * speed**2 / length**3,
        force / (rho * length**2),
    ]
    assert all(v == D(length=2, t=-2) for v in translation)
    shell = [sigma / length, chi / length, mus * epsilon * speed / length**2, mu * speed / length]
    assert all(v == p for v in shell)
    vdot = length**3 / time
    coefficient = rho * vdot**2
    assert coefficient == force * length**2
    assert coefficient / length**2 == force
    assert coefficient / (force * length**2) == D()  # Fig3 axis units divide dimensionally
    assert length**3 * (p / length) == force  # primary force
