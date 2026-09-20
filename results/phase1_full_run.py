import importlib.util
import sys
import time

LOG = open("results/phase1_full_run.txt", "w", encoding="utf-8", buffering=1)


def p(s):
    print(s)
    print(s, file=LOG)
    sys.stdout.flush()
    LOG.flush()


spec = importlib.util.spec_from_file_location(
    "ph1", "research/definitions/phase1_definitions_suite.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
if hasattr(m, "sys"):
    m.sys = sys

t0 = time.time()
m._results.clear()
m.REPORT = LOG
try:
    m.main()
except SystemExit as e:
    p(f"main() exited with code {e.code}")
p(f"wall_clock_total = {time.time() - t0:.2f} s")
LOG.close()