import freud
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# Параметры для расчета RDF
num_bins = 100  # Количество интервалов (бинов)
r_max = 8.0     # Максимальное расстояние для расчета RDF
# Инициализация объекта RDF
rdf = freud.density.RDF(bins=num_bins, r_max=r_max)

# Проход по всем кадрам
with gsd.hoomd.open('lip_mem/19_5/output.gsd', 'r') as f:
    for frame in f:
        type_ids = frame.particles.typeid  # Типы частиц
        mask = type_ids == 0  # Выбор частиц первого типа
        positions = frame.particles.position[mask]
        # positions = frame.particles.position
        box = freud.box.Box(*frame.configuration.box[:3])
        rdf.compute(system=(box, positions), reset=False)

# Получение усредненных результатов
r = rdf.bin_centers
g_r = rdf.rdf

# Визуализация
plt.plot(r, g_r, label='Averaged RDF')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Averaged Radial Distribution Function')
plt.legend()
plt.grid(True)
plt.show()
