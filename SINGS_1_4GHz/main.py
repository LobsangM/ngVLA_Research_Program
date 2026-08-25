import glob
import shutil
import subprocess
import sys
import os

scripts = [
    "SRF_MODEL_CASA.py",
    "ngVLA_config_A.py",
    "ngVLA_config_B.py",
    "VLA_config_A.py",
    #"VLA_config_B.py",
    "collage.py",
]

base = os.path.dirname(os.path.abspath(__file__))


def clean_casa_artifacts(base):
    print("\nLimpiando logs y temporales de CASA...")

    patterns = ["casa-*.log", "*.last", "temp_*.image", "temp_*.im"]

    for pattern in patterns:
        for path in glob.glob(os.path.join(base, pattern)):
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
            print("Eliminado:", path)


# `python main.py --clean` (o execfile con sys.argv conteniendo "--clean")
# solo limpia logs/temporales de CASA sin correr el pipeline (útil tras una
# corrida interrumpida/colgada). Sin guard de __name__ a propósito: este
# script también se lanza con execfile() dentro del shell de CASA, donde
# __name__ no vale "__main__".
if "--clean" in sys.argv:
    clean_casa_artifacts(base)
else:
    for script in scripts:
        print(f"\n{'='*50}")
        print(f"Ejecutando {script}...")
        print(f"{'='*50}\n")
        subprocess.run([sys.executable, script], cwd=base, check=True)

    clean_casa_artifacts(base)

    print("\nTodas las simulaciones terminadas y se realizaron los collages.")
