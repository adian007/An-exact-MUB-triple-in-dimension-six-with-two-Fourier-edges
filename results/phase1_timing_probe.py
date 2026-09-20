import importlib.util
import sys
import time

LOG = open("results/phase1_timing_probe.txt", "w", encoding="utf-8")


def p(s):
    print(s)
    print(s, file=LOG)
    LOG.flush()
    sys.stdout.flush()


spec = importlib.util.spec_from_file_location(
    "ph1", "research/definitions/phase1_definitions_suite.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

for name in ["test_T0_recipes", "test_T1_identity_fourier", "test_T2_gauge",
             "test_T3_redundancy", "test_T4_anchors", "test_T5_published_pair",
             "test_T6_equivalence", "test_T7_field"]:
    t = time.time()
    try:
        getattr(m, name)()
        p(f"### {name}: {time.time() - t:.2f} s")
    except Exception as e:
        p(f"### {name}: EXCEPTION after {time.time() - t:.2f} s: "
          f"{type(e).__name__}: {e}")
LOG.close()