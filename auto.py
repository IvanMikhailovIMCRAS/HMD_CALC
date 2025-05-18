from grafting import grafting
import os
from export_to_hoomd import gsd_compile
from calculation import calc
from calculation import calc_graft


if __name__ == '__main__':
    n = 3
    g = 2
    grafting(num_grafting=50, n=n, g=g, path_to_files='origin_files')
    # equlibration with frozen membrane
    os.system("mv ANGLS_exp ANGLS")
    os.system("mv BONDS_exp BONDS")
    os.system("mv COORD_exp COORD")
    os.system("mv FIXED_exp FIXED")
    os.system("cp origin_files/FIELD_frozen FIELD")
    os.system("cp origin_files/CONTR CONTR")
    os.system("./DPDL.exe")
    
    os.system("rm INFOR")
    os.system("rm TRACK")
    os.system("rm VELOCF")
    os.system("rm *.ent")
    os.system("rm energy.txt")
    
    os.system("mv COORDF COORD")
    bead_types = {0: "H", 1: "T", 2: "W", 3: "P"}
    gsd_compile(
        bead_types=bead_types,
        path="",
        name_output=f"init_n_{n}_g_{g}.gsd",
        )
    calc_graft(path=f"", label=f"_n_{n}_g_{g}")