import importlib.util
import sys
import time

INJECT = "D:\\MUBs in 6-dimension\\research\\definitions\\phase1_definitions_suite.py"
spec = importlib.util.spec_from_file_location("ph1", INJECT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

t0 = time.time()
for _ in range(3):
    try:
        m.test_T2_gauge()
    except Exception as e:
        break
    time.sleep(0.05)
print(f"### T2_gauge x3 wall={time.time()-t0:.2f}s")
