!######################################################################################
program DPD ! The program for Dissipative Particle Dynamic simulation 
!######################################################################################
	use CommonParam
	use ErrorList
	implicit none
	integer(4) N, NB, num_types ! число частиц, число связей, число типов частиц
	real(4) BOX(3) ! размер коробки
	real(4), allocatable, dimension(:) :: x, y, z, vx, vy, vz
	real(4), allocatable, dimension(:,:) :: FF ! "force-field"
	integer(4), allocatable, dimension(:) :: b1, b2, typ ! список связей, типы частиц
	integer(4), allocatable, dimension(:) :: list_fix    ! список фиксированных частиц
	integer(4), allocatable, dimension(:,:) :: list_ang  ! список "углов"
	integer(4) i, ioer, tmp_int, num_fix, num_ang
	character(1) tmp_char
	character(255) :: tmp_string
		
	! создаём файл для отчёта о выполнении кода и критических ошибках
	open(n_infor,file=name_infor)
	! проверяем существование в папке файла траектории TRACK, чтобы он не затёрся
	open(n_track,file=name_track,iostat=ioer,status='new')
	if (ioer.ne.0) call ERROR(1)
	close(n_track)
	! открываем файл с координатами (он должен существовать в папке)
	open(n_coord,file=name_coord,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(2)
	! считываем число частиц и размеры бокса
	read(n_coord,'(a)',advance='NO',iostat=ioer) tmp_string
	read(tmp_string,*,iostat=ioer) tmp_char, N, tmp_char, BOX(1), BOX(2), BOX(3)
	if (ioer.ne.0) then
		read(tmp_string,*,iostat=ioer) tmp_char, N, tmp_char, BOX(1)
		BOX(2) = BOX(1); BOX(3) = BOX(1)              
	endif
	if (ioer.ne.0.or.N.lt.3.or.BOX(1).le.0.0) call ERROR(3)
	write(n_infor,'(a,I0,a,F16.8)') 'N = ', N, '  BOX = ', BOX(1)
	! считываем координаты частиц
	allocate(x(1:N)); allocate(y(1:N)); allocate(z(1:N)) ! выделяем память под массивы
	allocate(typ(1:N))
	do i = 1, N
		read(n_coord,*,iostat = ioer) tmp_int, x(i), y(i), z(i), typ(i)
		if (ioer.ne.0 &
		& .or.abs(x(i)).gt.0.5*BOX(1).or.abs(y(i)).gt.0.5*BOX(2).or.abs(z(i)).gt.0.5*BOX(3)) then 
			write(n_infor,'(a,I0)') 'COORD file: string: ', i+1
			call ERROR(4)
		endif
	enddo
	close(n_coord)
	write(n_infor,'(a)') 'Coordinates have been read successfully.' 
	! открываем файл со скоростями (если его нет, задаём их все равными нулю) 
	open(n_veloc,file=name_veloc,iostat=ioer,status='old')
	allocate(vx(1:N)); allocate(vy(1:N)); allocate(vz(1:N)) ! выделяем память
	if (ioer.ne.0) then
		vx(:) = 0.0; vy(:) = 0.0; vz(:) = 0.0
		write(n_infor,'(a)') 'Warning! File VELOC does not exist! Velocities are set equal to zero.'
		write(n_infor,'(a)') 'Velocities are set equal to zero.'
	else
		read(n_veloc,*,iostat=ioer) tmp_char
		if (ioer.ne.0) call ERROR(5)
		do i = 1, N
			read(n_veloc,*,iostat = ioer) tmp_int, vx(i), vy(i), vz(i)
			if (ioer.ne.0) then 
				write(n_infor,'(a,I0)') 'VELOC file: string: ', i+1
				call ERROR(4)
			endif
		enddo
		close(n_veloc)
		write(n_infor,'(a)') 'Velocities have been read successfully.' 
	endif
	
	! открываем файл со списком связанных частиц (должен существовать)
	open(n_bonds,file=name_bonds,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(6)
	! считываем число связей
	read(n_bonds,*,iostat=ioer) tmp_char, NB
	if (ioer.ne.0.or.NB.lt.0) call ERROR(7)
	write(n_infor,'(a,I0)') 'NB = ', NB
	allocate(b1(1:NB)); allocate(b2(1:NB)) ! выделяем память под список связанных частиц
	do i = 1, NB
		read(n_bonds,*,iostat=ioer) b1(i), b2(i)
		if (ioer.ne.0.or.b1(i).ge.b2(i) &
		& .or.b1(i).lt.1.or.b1(i).gt.N.or.b2(i).lt.1.or.b2(i).gt.N) then 
			write(n_infor,'(a,I0)') 'BONDS file: string: ', i+1
			call ERROR(4)
		endif
	enddo
	close(n_bonds)
	write(n_infor,'(a)') 'Bonds have been read successfully.' 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_field,file=name_field,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(10)
	! считываем число связей
	read(n_field,*,iostat=ioer) tmp_char, num_types
	if (ioer.ne.0.or.num_types.lt.1) call ERROR(11)
	write(n_infor,'(a,I0)') 'num_types = ', num_types
	allocate(FF(1:num_types,1:num_types)) ! выделяем память под "силовое поле"
	read(n_field,*,iostat=ioer) FF
	if (ioer.ne.0) call ERROR(11)
	close(n_field)
	write(n_infor,'(a)') 'Force-field have been read successfully.' 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_fixed,file=name_fixed,iostat=ioer,status='old')
	if (ioer.ne.0) then
		num_fix = 0
	else
		! считываем число фиксированных атомов
		read(n_fixed,*,iostat=ioer) tmp_char, num_fix
		if (ioer.ne.0.or.num_fix.lt.1) call ERROR(12)
		write(n_infor,'(a,I0)') 'num_fixed = ', num_fix
		allocate(list_fix(1:N)); list_fix(:) = 0
		do i = 1, num_fix
			read(n_fixed,*,iostat=ioer) list_fix(i)
			if (ioer.ne.0.or.list_fix(i).gt.N.or.list_fix(i).lt.1) call ERROR(12)
		enddo
	endif	
	close(n_fixed)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_angls,file=name_angls,iostat=ioer,status='old')
	if (ioer.ne.0) then
		num_ang = 0
	else
		! считываем номера атомов включённых в валентный угол
		read(n_angls,*,iostat=ioer) tmp_char, num_ang
		if (ioer.ne.0.or.num_ang.lt.1) call ERROR(13)
		write(n_infor,'(a,I0)') 'num_ang = ', num_ang
		allocate(list_ang(1:3,1:N)); list_ang(:,:) = 0
		do i = 1, num_ang
			read(n_angls,*,iostat=ioer) list_ang(:,i)
			if (ioer.ne.0) call ERROR(13)
		enddo
	endif
	close(n_angls)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	write(n_infor,'(a)') 'Fixed points have been read successfully.' 
	
	call control(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, num_types, FF, &
	&            num_fix, list_fix, num_ang, list_ang)
	
	deallocate(x); deallocate(y); deallocate(z)
	deallocate(vx); deallocate(vy); deallocate(vz)
	deallocate(typ); deallocate(b1); deallocate(b2)
	deallocate(FF)
	if (num_fix.gt.0) deallocate(list_fix)
	write(n_infor,'(a)') 'Awesome! That worked out!' 
	close(n_infor)
	
	stop
end program DPD
!######################################################################################

!**************************************************************************************
subroutine control(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, num_types, FF, &
&                  num_fix, list_fix, num_ang, list_ang)
!**************************************************************************************
	use CommonParam
	use ErrorList
	use Engine
	use Lib
	implicit none
	integer(4), intent(in) :: N, NB, num_types, num_fix, num_ang
	real(4), intent(inout) :: BOX(3)
	real(4), intent(inout) :: x(1:N), y(1:N), z(1:N), vx(1:N), vy(1:N), vz(1:N)
	real(4), intent(inout) :: FF(1:num_types,1:num_types)
	integer(4), intent(inout) :: typ(1:N), list_fix(1:N), list_ang(1:3,1:num_ang)
	integer(4), intent(in) :: b1(1:NB), b2(1:NB)
	integer(4) num_step, num_snapshot, num_track, num_statis
	real(4) dt, l_bond, k_bond, k_ang
	integer(4) i, ioer
	integer(4) bond_list(1:N)                                             
	integer(4) nx, ny, nz, num_cell, N_dblist, N_dbcell 
	integer(4), allocatable, dimension(:) :: list1, list2, main_list, att_list
	
	open(n_contr,file=name_contr,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(8)
	
	read(n_contr,*,iostat=ioer) dt ! временной шаг интегрирования
	if (ioer.ne.0.or.dt.le.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) l_bond ! длина связи
	if (ioer.ne.0.or.l_bond.le.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) k_bond ! жёсткость связи
	if (ioer.ne.0.or.k_bond.lt.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) k_ang ! жёсткость угла
	if (ioer.ne.0.or.k_ang.lt.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer)	num_step ! число шагов
	if (ioer.ne.0.or.num_step.lt.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) num_snapshot ! через какое число шагов делать снимок
	if (ioer.ne.0.or.num_snapshot.lt.1) call ERROR(9)
	read(n_contr,*,iostat=ioer) num_track ! через какое число шагов печатать траекторию
	if (ioer.ne.0.or.num_track.lt.1) call ERROR(9)
	read(n_contr,*,iostat=ioer) delta_box ! с каким шагом растягиваем коробку по X и Y
	if (ioer.ne.0.or.delta_box.gt.(sqrt(BOX(1)*BOX(2)*BOX(3))-min(BOX(1),BOX(2)))/num_step) call ERROR(9)
		    	
	close(n_contr)
	write(n_infor,'(a)') 'CONTR-file have been read successfully.' 

	bond_list(:) = 0
	do i = 1, NB
		bond_list(b2(i)) = b1(i)     
	enddo

	nx = int(BOX(1)); ny = int(BOX(2)); nz = int(BOX(3))
	num_cell = nx * ny * nz   ! число ячеек
	N_dbcell = 13 * num_cell  ! число пар ячеек
	N_dblist = N * 256 ! число пар взаимодействующих атомов
	
	allocate(list1(0:N_dbcell));   allocate(list2(1:N_dbcell)) 
    allocate(main_list(num_cell)); allocate(att_list(N))
    list1(:) = 0; list2(:) = 0; main_list(:) = 0; att_list(:) = 0
       
    call DB_LIST(N_dbcell,BOX,list1,list2) ! находим список контактирующих ячеек
	
	call main(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, bond_list, dt, FF, &
	&        l_bond, k_bond, num_step, num_snapshot, num_track, N_dblist, num_fix, list_fix, &
	&        num_types, num_cell, N_dbcell, main_list, att_list, list1, list2, nx,ny,nz, &
	&        num_ang, list_ang, k_ang)
	
	deallocate(list1); deallocate(list2); deallocate(main_list); deallocate(att_list)
	
	return
end subroutine control
!**************************************************************************************
