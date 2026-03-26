"""Run a single test case for TP53."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_runner import run_test

result = run_test("TP53", "6M0J", 1)
print("\nDone.")
