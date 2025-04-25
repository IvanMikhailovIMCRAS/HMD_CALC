from export_to_hoomd import gsd_compile
from run_graft import calc_graft
from run import calc
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    bead_types = {0: "H", 1: "T", 2: "W", 3: "P"}
    for folder in ['g0', 'g3']:
        gsd_compile(bead_types=bead_types, path=f"new/{folder}", name_output=f"new/{folder}/init.gsd")
        calc_graft(f"new/{folder}")
        
    # bead_types = {0: "H", 1: "T", 2: "W"}
    # for folder in ['lip']:
    #     gsd_compile(bead_types=bead_types, path=f"new/{folder}", name_output=f"new/{folder}/init.gsd")
    #     calc(f"new/{folder}", free_steps=10000, calc_steps=30000)
    