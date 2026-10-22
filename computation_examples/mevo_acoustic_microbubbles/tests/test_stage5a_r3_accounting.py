"""Instrumentation must preserve values/exceptions and avoid timing double counts."""

from time import perf_counter
import pytest
from scripts.stage5a_r3.dynamic_audit import Accounting


def test_nested_partition_and_exception_unwinding():
    a = Accounting()
    start = perf_counter()
    with a.scope("outer"):
        with pytest.raises(ValueError):
            with a.scope("inner"):
                sum(range(1000))
                raise ValueError("candidate rejected")
        with a.scope("next"):
            sum(range(1000))
    wall = perf_counter() - start
    assert not a.stack
    assert dict(a.calls) == dict(outer=1, inner=1, next=1)
    assert all(v >= 0 for v in a.exclusive.values())
    assert sum(a.exclusive.values()) <= wall
    assert wall - sum(a.exclusive.values()) < 0.01


def test_qp_wrapper_preserves_result_and_failure():
    class Fake:
        def candidate(self, value):
            if value is None:
                raise ValueError("same rejected input")
            return value, {"unchanged": True}

    a = Accounting()
    table = []
    original = a.bind(Fake, "candidate", "native", table)
    try:
        token = object()
        result = Fake().candidate(token)
        assert result[0] is token
        with pytest.raises(ValueError, match="same rejected input"):
            Fake().candidate(None)
        assert table[0]["diagnostics"] == {"unchanged": True}
        assert table[1]["diagnostics"] is None
        assert table[1]["error"] is not None
        assert a.calls["native"] == 2
        assert not a.stack
    finally:
        Fake.candidate = original
