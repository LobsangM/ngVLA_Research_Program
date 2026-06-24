import subprocess
import sys
import os

scripts = [
    "SRF_MODEL_CASA.py",
    "ngVLA_config_A.py",
    "ngVLA_config_B.py",
    "VLA_config_A.py",
    "VLA_config_B.py",
]

base = os.path.dirname(os.path.abspath(__file__))

for script in scripts:
    print(f"\n{'='*50}")
    print(f"Ejecutando {script}...")
    print(f"{'='*50}\n")
    subprocess.run([sys.executable, script], cwd=base, check=True)

print("\nTodas las simulaciones terminadas.")
