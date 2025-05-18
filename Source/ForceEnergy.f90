!######################################################################################
module ForceEnergy
!######################################################################################
    use ErrorList
    use Lib
    use CommonParam
    use BoxMuller
    implicit none
		
Contains

!***************************************************************************************
subroutine calc_acceleration(N, NB, x, y, z, vx, vy, vz, BOX, FF, l_bond, &
&   k_bond, b1, b2, N_dblist, ub1, ub2, ax, ay, az, typ, num_types, sqrdt, num_ang, list_ang, k_ang, energy)
!***************************************************************************************
	integer(4), intent(in) :: N, NB, N_dblist, num_types, num_ang
	integer(4), intent(in) :: typ(1:N), list_ang(1:3,1:num_ang)
	real(4), intent(in) :: FF(1:num_types,1:num_types)
	integer(4), intent(in) :: b1(1:NB), b2(1:NB), ub1(0:N_dblist), ub2(1:N_dblist)
	real(4), intent(in) :: x(1:N), y(1:N), z(1:N), vx(1:N), vy(1:N), vz(1:N)
	real(4), intent(in) :: sqrdt
	real(4), intent(in) :: l_bond, k_bond, k_ang, BOX(3)
	real(4), intent(out) :: ax(1:N), ay(1:N), az(1:N), energy
	integer(4) i
	real(4) dx, dy, dz, dvx, dvy, dvz, r, force, force_C, force_R, force_D, omega, ksi
	real(4) f_ang(1:3,1:3), gm, th, Cs
	real(4) kBT
	
	kBT = 1.0
	energy = 0.0

	ax(:) = 0.0; ay(:) = 0.0; az(:) = 0.0
	f_ang(:,:) = 0.0
	
        ! рассчитываем НЕвалентные взаимодействия
	do i = 1, ub1(0)
		dx = x(ub1(i))-x(ub2(i)); dy = y(ub1(i))-y(ub2(i)); dz = z(ub1(i))-z(ub2(i))
		call periodic_test(dx,BOX(1));	call periodic_test(dy,BOX(2));	call periodic_test(dz,BOX(3))
		r = sqrt(dx**2 + dy**2 + dz**2)
		if (r.gt.r_cut) then
			force = 0.0
		else
			omega = 1.0 - r/r_cut
			call GaussRand(ksi)
			dvx = vx(ub1(i))-vx(ub2(i)); dvy = vy(ub1(i))-vy(ub2(i)); dvz = vz(ub1(i))-vz(ub2(i))	
			
			if (FF(typ(ub1(i)),typ(ub2(i))).gt.34.999) then
				if (FF(typ(ub1(i)),typ(ub2(i))).lt.74.999) then
				    gm = 9.0
					th = 4.24264069
				else
					gm = 20.0
					th = 6.32455532
				endif
			else
				gm = gamma
				th = theta
			endif
					
			force_C = FF(typ(ub1(i)),typ(ub2(i))) * omega / r
			energy = energy + 0.5*FF(typ(ub1(i)),typ(ub2(i)))*(1.0- r)**2
			force_R = th * omega * ksi / sqrdt / r
			force_D = - gm * omega * omega * (dx*dvx + dy*dvy + dz*dvz) / r / r	
					
			force = force_C + force_R + force_D
			force = force * kBT
		endif
		ax(ub1(i)) = ax(ub1(i)) + force * dx; ax(ub2(i)) = ax(ub2(i)) - force * dx 
		ay(ub1(i)) = ay(ub1(i)) + force * dy; ay(ub2(i)) = ay(ub2(i)) - force * dy
		az(ub1(i)) = az(ub1(i)) + force * dz; az(ub2(i)) = az(ub2(i)) - force * dz
	enddo

	! рассчитываем валентные взаимодействия
	
	do i = 1, NB
		dx = x(b1(i))-x(b2(i)); dy = y(b1(i))-y(b2(i)); dz = z(b1(i))-z(b2(i))
		call periodic_test(dx,BOX(1));	call periodic_test(dy,BOX(2));	call periodic_test(dz,BOX(3))
		r = sqrt(dx**2 + dy**2 + dz**2)
		if (typ(b2(i)).eq.4) then
			force_C = 4222.6 * (0.4125/r - 1.0) 
		else
			force_C = k_bond * (l_bond/r - 1.0) 
		endif
		
		
		if (r.gt.r_cut) then
			force = force_C
		else
			omega = 1.0 - r/r_cut
			call GaussRand(ksi)
			dvx = vx(b1(i))-vx(b2(i)); dvy = vy(b1(i))-vy(b2(i)); dvz = vz(b1(i))-vz(b2(i))	
			
			if (FF(typ(b1(i)),typ(b2(i))).gt.34.999) then
				if (FF(typ(b1(i)),typ(b2(i))).lt.74.999) then
				    gm = 9.0
					th = 4.24264069
				else
					gm = 20.0
					th = 6.32455532
				endif
			else
				gm = gamma
				th = theta
			endif

			force_R = th * omega * ksi / sqrdt / r
			force_D = - gm * omega * omega * (dx*dvx + dy*dvy + dz*dvz) / r / r	
					
			force = force_C + force_R + force_D
			force = force * kBT
		endif
		ax(b1(i)) = ax(b1(i)) + force * dx; ax(b2(i)) = ax(b2(i)) - force * dx 
		ay(b1(i)) = ay(b1(i)) + force * dy; ay(b2(i)) = ay(b2(i)) - force * dy
		az(b1(i)) = az(b1(i)) + force * dz; az(b2(i)) = az(b2(i)) - force * dz
	enddo
	
	! валентные углы
	do i = 1, num_ang
		call angular_projections(x(list_ang(1,i)), x(list_ang(2,i)), x(list_ang(3,i)), &
&                                y(list_ang(1,i)), y(list_ang(2,i)), y(list_ang(3,i)), &
&                                z(list_ang(1,i)), z(list_ang(2,i)), z(list_ang(3,i)), &
&                                f_ang(1,1),f_ang(1,2),f_ang(1,3), &
&                                f_ang(2,1),f_ang(2,2),f_ang(2,3), &
&                                f_ang(3,1),f_ang(3,2),f_ang(3,3), Cs) 
		f_ang(:,:) = f_ang(:,:) * kBT
		if (typ(list_ang(2,i)).eq.4) then
			ax(list_ang(1,i)) = ax(list_ang(1,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(1,1) 
			ax(list_ang(2,i)) = ax(list_ang(2,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(1,2)
			ax(list_ang(3,i)) = ax(list_ang(3,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(1,3)
			ay(list_ang(1,i)) = ay(list_ang(1,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(2,1)
			ay(list_ang(2,i)) = ay(list_ang(2,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(2,2)
			ay(list_ang(3,i)) = ay(list_ang(3,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(2,3)
			az(list_ang(1,i)) = az(list_ang(1,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(3,1)
			az(list_ang(2,i)) = az(list_ang(2,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(3,2)
			az(list_ang(3,i)) = az(list_ang(3,i)) - 32.9892 * (Cs - 0.6427876) * f_ang(3,3) 
		else
			ax(list_ang(1,i)) = ax(list_ang(1,i)) + k_ang * f_ang(1,1) 
			ax(list_ang(2,i)) = ax(list_ang(2,i)) + k_ang * f_ang(1,2)
			ax(list_ang(3,i)) = ax(list_ang(3,i)) + k_ang * f_ang(1,3)
			ay(list_ang(1,i)) = ay(list_ang(1,i)) + k_ang * f_ang(2,1)
			ay(list_ang(2,i)) = ay(list_ang(2,i)) + k_ang * f_ang(2,2)
			ay(list_ang(3,i)) = ay(list_ang(3,i)) + k_ang * f_ang(2,3)
			az(list_ang(1,i)) = az(list_ang(1,i)) + k_ang * f_ang(3,1)
			az(list_ang(2,i)) = az(list_ang(2,i)) + k_ang * f_ang(3,2)
			az(list_ang(3,i)) = az(list_ang(3,i)) + k_ang * f_ang(3,3)
		endif
	enddo
	
return
end subroutine calc_acceleration
!***************************************************************************************
    
!######################################################################################
end module ForceEnergy
!######################################################################################
