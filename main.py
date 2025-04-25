from export_to_hoomd import gsd_compile
from run import calc
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # bead_types = {0: "H", 1: "T", 2: "W"}
    # for folder in ['19_5']:
    #     gsd_compile(bead_types=bead_types, path=f"lip_mem/{folder}", name_output=f"lip_mem/{folder}/init.gsd")
    #     calc(f"lip_mem/{folder}")
    gamma = []
    gamma_err = []
    for folder in ['17', '18', '19', '19_5', '20', '21', '21_5']:
        tension = np.loadtxt(f"lip_mem/{folder}/tension.txt")
        gamma.append(np.mean(tension))
        # gamma_err.append(abs(np.mean(tension[:len(tension)//2]) - np.mean(tension[len(tension)//2:]))/2)
        gamma_err.append(np.sqrt(np.std(gamma)))
        
    # Построение графика с интервалами ошибок
    L_box = np.array([17, 18, 19, 19.5, 20, 21, 21.5])
    plt.errorbar(L_box**2/530*2, gamma, yerr=gamma_err, fmt='-o', capsize=5)

    # Настройка графика
    plt.xlabel('$a_{prj}~(r_0^2)$', fontsize=14)
    plt.ylabel('$\Sigma~(k_B T/r_0^2)$', fontsize=14)
    plt.title('Поверхностное натяжение "голой" липидной мембраны')
    plt.grid(True)

    # Отображение графика
    plt.savefig("tention.jpg", dpi=300)