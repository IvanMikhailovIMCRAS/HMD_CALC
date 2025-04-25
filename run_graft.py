import hoomd
from hoomd import md
import numpy as np
import time
import os

def calc_graft(path: str, free_steps: int=100000, calc_steps: int=1e6, trigger: int=1000) -> None:
    # Параметры потенциала угла
    k = 15.0  # Жесткость потенциала
    num_points = 100  # Количество точек в таблице
    theta_min = 0.0   # Минимальный угол (в радианах)
    theta_max = np.pi # Максимальный угол (в радианах)

    # Создание массивов для ВАЛЕНТНОГО угла, энергии и силы
    theta = np.linspace(theta_min, theta_max, num_points)
    energy = k * (1 + np.cos(theta))  # Энергия: V(theta) = k * (1 + cos(theta))
    force = k * np.sin(theta)        # Сила: F(theta) = -dV/dtheta = k * sin(theta)
    
    k2 = 2*16.4946
    theta0 = 130.0
    energy2 = k2 * (np.cos(theta0) - np.cos(theta))**2
    force2 = - 2.0 * k2 * (np.cos(theta0) - np.cos(theta)) * np.sin(theta)

    # Добавление силы через таблицу
    table_potential = hoomd.md.angle.Table(width=num_points)
    table_potential.params['HTT'] = dict(U=energy, tau=force)
    table_potential.params['TTT'] = dict(U=energy, tau=force)
    table_potential.params['HPP'] = dict(U=energy2, tau=force2)
    table_potential.params['PPP'] = dict(U=energy2, tau=force2)
    
    # Буфер для списка взаимодействующих частиц (40% от общего числа)
    nl = hoomd.md.nlist.Cell(buffer = 0.4)
    
    # Параметры DPD взаимодействия между несвязанными бидами
    dpd = md.pair.DPD(default_r_cut=1.0, kT=1.0, nlist=nl)
    dpd.params[('H', 'H')] = dict(A=30.0, gamma = 4.5)
    dpd.params[('H', 'T')] = dict(A=35.0, gamma = 9.0)
    dpd.params[('H', 'W')] = dict(A=30.0, gamma = 4.5)
    dpd.params[('T', 'T')] = dict(A=10.0, gamma = 4.5)
    dpd.params[('T', 'W')] = dict(A=75.0, gamma = 20.0)
    dpd.params[('W', 'W')] = dict(A=25.0, gamma = 4.5)
    
    dpd.params[('P', 'P')] = dict(A=25.0, gamma = 4.5)
    dpd.params[('P', 'H')] = dict(A=26.3, gamma = 4.5)
    dpd.params[('P', 'T')] = dict(A=33.7, gamma = 4.5)
    dpd.params[('P', 'W')] = dict(A=26.3, gamma = 4.5)      

    # Гармонический потенциал связей в липидах
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['HH'] = dict(k=128.0, r0=0.5)
    harmonic.params['HT'] = dict(k=128.0, r0=0.5)
    harmonic.params['TT'] = dict(k=128.0, r0=0.5)
    
    harmonic.params['HP'] = dict(k=2*2111.3, r0=0.4125)
    harmonic.params['PP'] = dict(k=2*2111.3, r0=0.4125)

    # Подключаем GPU девайс (для распараллеливания на NVIDIA)
    gpu = hoomd.device.GPU()
    sim = hoomd.Simulation(device=gpu, seed=200)
    
    # Считываем начальную конфигурацию системы
    sim.create_state_from_gsd(filename=os.path.join(path, "init.gsd"))
    # При использовании gpu без этой строчки не вычисляется тензор давления!!!
    sim.always_compute_pressure = True
    # шаг интегрирования (TODO half step hook)
    integrator = hoomd.md.Integrator(dt=0.01)
    # выбираем NVT ансамбль при постоянном объёме
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
    sim.operations.integrator = integrator
    sim.operations.integrator.methods = [nvt]
    # добавляем все типы сил в итератор
    integrator.forces.append(dpd)
    integrator.forces.append(harmonic)
    integrator.forces.append(table_potential)
    # integrator.forces.append(angl_potential)

    # Запрашиваем нужные парамтры системы 
    box = sim.state.box  # Получаем объект Box
    Lx, Ly, Lz = box.Lx, box.Ly, box.Lz  # Размеры ящика
    V = Lx * Ly * Lz

    # сюда будем сохранять поверхностное натяжение
    gamma = []

    # засекаем секундомер
    t1 = time.time()
    
    # прогоняем free_steps для уравновешивания системы
    sim.run(free_steps)

    # Настройка вывода данных
    gsd = hoomd.write.GSD(trigger=trigger, filename=os.path.join(path,'output.gsd'))
    sim.operations.writers.append(gsd)
    print(calc_steps//trigger)
    for _ in range(int(calc_steps)//trigger):
        sim.run(trigger)
        # thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
        thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.Type(['H','T']))
        sim.operations.computes.append(thermo)
        p_xx = thermo.pressure_tensor[0]
        p_yy = thermo.pressure_tensor[3]
        p_zz = thermo.pressure_tensor[5]
        gamma.append((p_zz - 0.5*(p_xx+p_yy))*V/Lx/Ly)
        
    # останавливаем секундомер
    t2 = time.time()
    print(f"lead time: {t2-t1} sec")
    
    # печатаем значения поверхностного натяжения в файл
    tension = np.array(gamma)
    np.savetxt(os.path.join(path, 'tension.txt'), tension)
    