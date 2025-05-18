!######################################################################################
module Output
!######################################################################################
    use CommonParam
    implicit none
		
Contains

!***************************************************************************************
subroutine show(N, NB, b1, b2, BOX, x, y, z, typ, num) ! печатает ent-файл
!***************************************************************************************
	integer(4),intent(in) :: N, NB, num
	integer(4), intent(in) :: b1(1:NB), b2(1:NB), typ(1:N)
	real(4), intent(in) :: BOX(3)
	real(4), intent(in) :: x(1:N), y(1:N), z(1:N)
	character(128) Name_file
	real(4) dx, dy, dz
	integer(4) i
	
		
    ! создаём имя файла Name_file в формате "picture_000000000123.ent"
	write(Name_file,'(a,I0.12,a)')'picture_',num,'.ent'
	
	open(n_ent,file = Name_file)
	! координаты шариков:
	do i = 1, N
		if (typ(i).gt.6) then
			write(n_ent,'(a,I5,2x,a,I12,4x,F8.3,F8.3,F8.3)') &
			&           'HETATM',i,palette(6:6),i,x(i),y(i),z(i)
		else
			if (typ(i).eq.3) then
				write(n_ent,'(a,I5,2x,a,I12,4x,F8.3,F8.3,F8.3)') &
				&           'HETATM',0,palette(typ(i):typ(i)),0,x(i),y(i),z(i)
			else
				write(n_ent,'(a,I5,2x,a,I12,4x,F8.3,F8.3,F8.3)') &
				&           'HETATM',i,palette(typ(i):typ(i)),i,x(i),y(i),z(i)
			endif
		endif
	enddo
	! список связей
	do i = 1, NB
		dx = x(b1(i))-x(b2(i)); dy = y(b1(i))-y(b2(i)); dz = z(b1(i))-z(b2(i))
		if (abs(dx).lt.BOX(1)*0.5.and.abs(dy).lt.BOX(2)*0.5.and.abs(dz).lt.BOX(3)*0.5) then
			write(n_ent,'(a,I5,I5)') 'CONECT', b1(i), b2(i)
		endif
	enddo
	
	close(n_ent)

return
end subroutine show
!***************************************************************************************

!**************************************************************************************
subroutine print_track(N, x, y, z, typ, num_step, BOX)
!**************************************************************************************
implicit none
    integer(4), intent(in) :: N, num_step
    real(4), intent(in) :: x(1:N), y(1:N), z(1:N)
    integer(4), intent(in) :: typ(1:N) 
    real(4), intent(in) :: BOX(3)
    integer(4) i
    logical store_label_track 
    save store_label_track ! делаем метку о первом вызове процедуры статической
    data store_label_track / .true. /
    
    open(n_track,file=name_track,access='append',status='old')
    if (store_label_track) then
    	write(n_track,'(a,1x,I0,1x,a,1x,3(E16.8))') 'num_atoms', N, 'box_size', BOX
    	store_label_track = .false.
    endif
    
    write(n_track,'(a,1x,I0,1x,a,1x,3(E16.8))') 'step', num_step, 'box_size', BOX
    do i = 1, N
    	write(n_track,'(I12,3(1xE16.8),1x,I0)') i, x(i), y(i), z(i), typ(i)
    enddo
    
    close(n_track)    
return
end subroutine print_track
!**************************************************************************************

!**************************************************************************************
subroutine print_coordf(N, x, y, z, typ, BOX)
!**************************************************************************************
implicit none
    integer(4), intent(in) :: N
    real(4), intent(in) :: x(1:N), y(1:N), z(1:N)
    integer(4), intent(in) :: typ(1:N)
    real(4), intent(in) :: BOX(3)
    integer(4) i
    
    open(n_coordf,file=name_coordf)
    
	write(n_coordf,'(a,1x,I0,1x,a,1x,3(E16.8))') 'num_atoms', N, 'box_size', BOX
	do i = 1, N
		write(n_coordf,'(I12,3(1xE16.8),1x,I0)') i, x(i), y(i), z(i), typ(i)
	enddo
	
	close(n_coordf)
      
return
end subroutine print_coordf
!**************************************************************************************

!**************************************************************************************
subroutine print_velocf(N, vx, vy, vz)
!**************************************************************************************
implicit none
    integer(4), intent(in) :: N
    real(4), intent(in) :: vx(1:N), vy(1:N), vz(1:N)
    integer(4) i
    
    open(n_velocf,file=name_velocf)
    
	write(n_velocf,'(a,1x,I0)') 'num_atoms', N
	do i = 1, N
		write(n_velocf,'(I12,3(1xE16.8))') i, vx(i), vy(i), vz(i)
	enddo
	
	close(n_velocf)
      
return
end subroutine print_velocf
!**************************************************************************************
    
!######################################################################################
end module Output
!######################################################################################
