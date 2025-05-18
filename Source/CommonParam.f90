!######################################################################################
module CommonParam 
!######################################################################################
	implicit none
!!!!!!!!! descriptors for input files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_coord = 11
	character(5), public, parameter :: name_coord = 'COORD' ! координаты частиц
	integer(4), public, parameter :: n_bonds = 12
	character(5), public, parameter :: name_bonds = 'BONDS' ! список связанных частиц
	integer(4), public, parameter :: n_veloc = 13
	character(5), public, parameter :: name_veloc = 'VELOC' ! начальные скорости
	integer(4), public, parameter :: n_contr = 14
	character(5), public, parameter :: name_contr = 'CONTR' ! управляющие параметры
	integer(4), public, parameter :: n_field = 15
	character(5), public, parameter :: name_field = 'FIELD' ! "силовое поле"
	integer(4), public, parameter :: n_fixed = 16
	character(5), public, parameter :: name_fixed = 'FIXED' ! фиксированные атомы
	integer(4), public, parameter :: n_angls = 17
	character(5), public, parameter :: name_angls = 'ANGLS' ! список "углов"	
!!!!!!!!! descriptors fot output files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_ent = 20 ! дескриптор для выходных ent-файлов
	integer(4), public, parameter :: n_track = 21
	character(5), public, parameter :: name_track = 'TRACK' ! печать траектории
	integer(4), public, parameter :: n_coordf = 22
	character(6), public, parameter :: name_coordf = 'COORDF' ! финальная конфигурация
	integer(4), public, parameter :: n_velocf = 23
	character(6), public, parameter :: name_velocf = 'VELOCF' ! финальные скорости
	integer(4), public, parameter :: n_infor = 24
	character(5), public, parameter :: name_infor = 'INFOR'   ! отчёт о выполнении кода
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	character(6), public, parameter :: palette = 'NOCSPF'
	real(4), public, parameter :: theta = 3.0
	real(4), public, parameter :: gamma = 4.5
	real(4), public, parameter :: r_cut = 1.0
	real(4), public, parameter :: lambda = 0.65
	
	real(4) box_x, box_y, box_volume, delta_box
	
Contains



end module CommonParam



