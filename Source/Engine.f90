!######################################################################################
module Engine
!######################################################################################
	use CommonParam
	use ErrorList
	use Output
	use ForceEnergy
	use Lib
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
Contains

!***************************************************************************************
subroutine main(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, bond_list, dt, FF, &
&        l_bond, k_bond, num_step, num_snapshot, num_track, N_dblist, num_fix, list_fix,  &
&        num_types, num_cell, N_dbcell, main_list, att_list, list1, list2, nx,ny,nz, &
&        num_ang, list_ang, k_ang)
!***************************************************************************************
	implicit none
	integer(4), intent(in) :: N, NB, num_types, nx,ny,nz, num_fix, num_ang
	real(4), intent(inout) :: BOX(3)
	real(4), intent(inout) :: x(1:N), y(1:N), z(1:N), vx(1:N), vy(1:N), vz(1:N)
	real(4), intent(inout) :: FF(1:num_types,1:num_types)
	integer(4), intent(inout) :: typ(1:N), list_fix(1:N), list_ang(1:3,1:num_ang)
	integer(4), intent(in) :: b1(1:NB), b2(1:NB)
	integer(4), intent(in) :: num_step, num_snapshot, num_track
	integer(4), intent(in) :: N_dblist, num_cell, N_dbcell
	real(4), intent(in) ::  dt, l_bond, k_bond, k_ang
	integer(4), intent(in) ::  bond_list(1:N), list1(0:N_dbcell), list2(1:N_dbcell)
	integer(4), intent(inout) :: main_list(1:num_cell), att_list(1:N)
	
	real(4) ax(1:N),ay(1:N),az(1:N)
	real(4) ax_old(1:N),ay_old(1:N),az_old(1:N)
	real(4) vx_old(1:N),vy_old(1:N),vz_old(1:N)
	real(4) xf(1:N),yf(1:N),zf(1:N)  
	integer(4) step, i
	integer(4) ub1(0:N_dblist),ub2(1:N_dblist)
	real(8) t(3)
	real(4) sqrdt, energy
	
	xf(:) = x(:); yf(:) = y(:); zf(:) = z(:)
	
	sqrdt = sqrt(dt)
	
	! delta_box = 0.01
	box_x = BOX(1)
	box_y = BOX(2)
	box_volume = BOX(1) * BOX(2) * BOX(3)
	
	ax(:) = 0.0; ay(:) = 0.0; az(:) = 0.0 
	ax_old(:) = 0.0; ay_old(:) = 0.0; az_old(:) = 0.0 
	vx_old(:) = 0.0; vy_old(:) = 0.0; vz_old(:) = 0.0 
	ub1(:) = 0; ub2(:) = 0
	
	open(33, file = 'energy.txt')
			
	call attached_list(num_cell, N, x, y, z, nx, ny, nz, BOX, main_list, att_list)	
	call scan_list(N, num_cell, N_dblist, N_dbcell, main_list, att_list, list1, list2,  &
&                    bond_list, ub1, ub2)
		
	! печатаем начальную конфигурацию
	call show(N, NB, b1, b2, BOX, x, y, z, typ, 0) 
			
	! рассчитываем начальные ускорения	
	call calc_acceleration(N, NB, x, y, z, vx, vy, vz, BOX, FF, l_bond, k_bond, &
&		b1, b2, N_dblist, ub1, ub2, ax, ay, az, typ, num_types, sqrdt, num_ang, list_ang, k_ang, energy) 

	! дважды запускаем секундомер
	call cpu_time(t(1)); call cpu_time(t(2))
	!++++++++++++++++++++++ НАЧАЛО ОСНОВНОГО ЦИКЛА ++++++++++++++++++++++++++++++++++++
	do step = 1, num_step
		BOX(1) = box_x + delta_box * step
		BOX(2) = box_y + delta_box * step
		BOX(3) = box_volume / BOX(1) / BOX(2)
		! рассчитываем новые координаты x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
		x(:) = x(:) + vx(:)*dt + 0.5*ax(:)*dt**2
		y(:) = y(:) + vy(:)*dt + 0.5*ay(:)*dt**2
		z(:) = z(:) + vz(:)*dt + 0.5*az(:)*dt**2
		! если атомы фиксированные - возвращаем всё на место
		do i = 1, num_fix
			x(list_fix(i)) = xf(list_fix(i)) 
			y(list_fix(i)) = yf(list_fix(i))
			z(list_fix(i)) = zf(list_fix(i)) 
		enddo
		! проверяем на выход за пределы коробки
		do i = 1, N
			call periodic_test(x(i),BOX(1))
			call periodic_test(y(i),BOX(2))
			call periodic_test(z(i),BOX(3))
		enddo
	
		call attached_list(num_cell, N, x, y, z, nx, ny, nz, BOX, main_list, att_list)	
		call scan_list(N, num_cell, N_dblist, N_dbcell, main_list, att_list, list1, list2,  &
&                     bond_list, ub1, ub2)

		vx_old(:) = vx(:)
		vy_old(:) = vy(:)
		vz_old(:) = vz(:)

		vx(:) = vx(:) + lambda*ax(:) * dt 
		vy(:) = vy(:) + lambda*ay(:) * dt 
		vz(:) = vz(:) + lambda*az(:) * dt 
	
		ax_old(:) = ax(:)
		ay_old(:) = ay(:)
		az_old(:) = az(:)
	
		call calc_acceleration(N, NB, x, y, z, vx, vy, vz, BOX, FF, l_bond, k_bond, &
&		b1, b2, N_dblist, ub1, ub2, ax, ay, az, typ, num_types, sqrdt, num_ang, list_ang, k_ang, energy) 
		

		vx(:) = vx_old(:) + 0.5*(ax(:) + ax_old(:)) * dt 
		vy(:) = vy_old(:) + 0.5*(ay(:) + ay_old(:)) * dt 
		vz(:) = vz_old(:) + 0.5*(az(:) + az_old(:)) * dt 
			
		!!!!!! ВЫВОДИМ НА ПЕЧАТЬ ВСЯКОЕ РАЗНОЕ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		! печатаем снимки системы каждые num_snapshot шагов
		if(modulo(step,num_snapshot).eq.0) call show(N,NB,b1,b2,BOX,x,y,z,typ,step)
		
		! печатаем траекторию каждые num_track шагов
		if(modulo(step,num_track).eq.0) then 
			call print_track(N, x, y, z, typ, step, BOX)	
			write(33, *) step, energy
		endif			
	enddo
	close(33)
	!++++++++++++++++++++++ КОНЕЦ ОСНОВНОГО ЦИКЛА +++++++++++++++++++++++++++++++++++++
	
	! запускаем в третий раз секундомер и подсчитываем время работы основного цикла
	call cpu_time(t(3))
	write(n_infor,*) 'Lead time (sec): ', t(3) - t(2) - (t(2) - t(1))
	
	call print_coordf(N, x, y, z, typ, BOX)  ! печатаем финальную конфигурацию
	call print_velocf(N, vx, vy, vz)         ! печатаем финальные скорости
	call show(N,NB,b1,b2,BOX,x,y,z,typ,num_step+1) ! фотокарточка на прощание

return
end subroutine main
!**************************************************************************************

!######################################################################################
end module Engine
!######################################################################################


