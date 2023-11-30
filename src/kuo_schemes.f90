program kuo_schemes
 implicit none
 integer, parameter :: nlev=60, nprof=1535, klev_max=20
 real, parameter :: Lv0=2.501E6, Cp=1005.0, Rd=287.06 
 real, parameter :: Eps=5.0E-5, Cpv=1846.0, Cl=4218.0, rg=9.81
! Values above are consistent with Marquet and Malardel + ECMWF constants 
! Other constants proposed by Emanuel (1994) : Cpv=1870 - Cl=4190 
 real, dimension(nlev) :: p, t, qv, gz, w, u, v, tc, qvc, ptent, ptenq, &
                       &  dtdt, dqdt, qwb, twb, rh, alpha_env, Cps_env, tve
 real ::  bkuo, Cps, cvgh, cvgu, dt, dummy, dz, gz2, p2, ps, qc, qv1, qv2, qw, rain, &
        & rhmean, t1, t2, ts, zlat, zlon, zlsm, tvc, zcape, zcv1, zcv2, zint, zint2, &
        & zint_t, zint_q, zlam, zmean, s1, tw, qsat, dqsat, alpha_w, Lh, w_lcl,  &
        & t_excess, q_excess 
 integer :: kdummy, i, irec, k                      
 real, dimension(nlev+1) :: ph, gzh
 integer, dimension(nlev) :: icond
 logical :: lsub
 character*3 :: expid
 
 expid='004'
!
 open (unit=20,file='../data_in/profiles.dat',status='old')
 open (unit=21,file='../data_in/input_model_mean_values.dat',status='old')
 open (unit=200,file='../data_out/conv_tendencies_kuo65_'//expid//'.dat')
 open (unit=201,file='../data_out/conv_tendencies_kuogeleyn_'//expid//'.dat')
 open (unit=202,file='../data_out/conv_tendencies_kuosymmetric_'//expid//'.dat')
 open (unit=203,file='../data_out/conv_tendencies_kuoanthes_'//expid//'.dat')
 open (unit=300,file='../data_out/atmospheric_profiles_'//expid//'.dat')
!
! Model time step
!
 dt = 900.0
!
! Logical to account for latent heat of sublimation 
! 
 lsub = .false.
!
! Conversion factors for heating and moistening rates
! 
 zcv1 = 86400.0
 zcv2 = zcv1*1.0E3
!
! Possibility of modifying (T,q) properties of the parcel (from surface fluxes)
! 
 t_excess = 0.0
 q_excess = 0.0
!
! Read atmospheric profiles
!  
 do irec=1,900
   read(20,*)
   read(20,*) zlat, zlon, zlsm, ts, ps
   read(20,*)
   do k=1,nlev
     read(20,*) t(k),qv(k),p(k),ph(k),w(k),u(k),v(k)
     rh(k) = qv(k)/qsat(p(k),t(k))
     Cps_env(k) = Cp*(1.0 - qv(k)) + qv(k)*Cpv
     alpha_env(k) = Lh(t(k),lsub)/Cps_env(k)
   enddo
   ph(nlev+1) = ps  
 enddo 
! 
 close(unit=20)
!
! Read dynamical tendencies for T and q
!
 do k=1,nlev
   read(21,*) kdummy,dummy,dummy,ptent(k),ptenq(k)
 enddo 
 ptent = 0.25*ptent 
! 
 close(unit=21)
!
! Compute geopotential (integrate hydrostatic equation)
!
 gzh(nlev+1) = 0.0
 do k=nlev,1,-1
   tve(k) = t(k)*(1. + 0.608*qv(k))
   gzh(k) = gzh(k+1) + Rd*tve(k)*(ph(k+1)- ph(k))/p(k)
   gz(k) = 0.5*(gzh(k+1) + gzh(k))
 enddo 
!
! Define moist pseudo adiabat (condensed water is removed after being created) 
!
 t1 = t(nlev) + t_excess
 qv1 = qv(nlev) + q_excess
 t2 = t1
 qv2 = qv1 
 qc = 0.0
 icond(:) = 0
 Cps = Cp*(1.0 - qv1) + qv1*Cpv
 s1 = Cps*t1 + gz(nlev) + Lh(t1,lsub)*qv1
!
! Initialise cloud values to environmental properties
! 
 tc = t
 qvc = qv 
 twb = t
 qwb = qv
! 
 zcape = 0.0
 do k=nlev-1,klev_max,-1
   dz = (gz(k+1) - gz(k))/rg
   p2 = p(k)
   gz2 = gz(k)
   t2 = (s1 - gz2 - Lh(t2,lsub)*qv1)/Cps
   t1 = t2 ! save this temperature for the "no saturation" case
   qv2 = qv1 
   tw = t(k)
   alpha_w = Lh(tw,lsub)/(Cp*(1.0 - qsat(p2,tw)) + qsat(p2,tw)*Cpv)
   zlam = (s1 - gz2) - Eps*dz*(Cps_env(k)*t(k) + Lh(t(k),lsub)*qv(k))
   do i=1,5
     t2 = t2 - ((Cps*t2 + Lh(t2,lsub)*qsat(p2,t2))*(1.0 - Eps*dz) - zlam)/ &
        & ((Cps + Lh(t2,lsub)*dqsat(p2,t2) + (Cpv - Cp)*dqsat(p2,t2)*t2 + & 
        & (Cpv - Cl)*qsat(p2,t2))*(1.0 - Eps*dz))           
     tw = tw - (tw + alpha_w*qsat(p2,tw) - t(k) - alpha_w*qv(k))/ &
        & (1.0 + alpha_w*dqsat(p2,tw))   
   enddo 
   qw = qsat(p2,tw)
   qwb(k) = qw
   twb(k) = tw
   if (qv2 > qsat(p2,t2) .and. qv2 > 1.E-9) then
     if (icond(k) == 0 .and. icond(k+1) == 0) w_lcl = w(k)
     icond(k) = 1
     qc  = qv2 - qsat(p2,t2)
     qv1 = qsat(p2,t2)
   else
     icond(k) = 0
     qc  = 0.0
     qv1 = qv2 
     t2 = t1
   endif 
!
!  Estimate regions with positive buoyancy from virtual temperature
!  Activate convection only when the vertical velocity is positive at LCL
!    
   tvc = t2*(1.0 + 0.608*qv1 - qc) ! inclusion of water loading in virtual temperature
   if (tvc > tve(k) .and. icond(k) == 1) then
     if (w_lcl > 0.0) icond(k) = 2  
     zcape = zcape - rg*dz*(tvc - tve(k))/tve(k)
   endif 
!
! Compute new moist static energy at level (k+1) for level (k)
!   
   Cps = Cp*(1.0 - qv1) + qv1*Cpv
   s1 = Cps*t2 + gz2 + Lh(t2,lsub)*qv1
!
!  Save T and q in-cloud properties
!   
   qvc(k) = qv1
   tc(k) = t2
 enddo 
!
! Moist static energy integral + Convergence of humidity and MSE
! 
 cvgu = 0.0
 cvgh = 0.0
 zint = 0.0
 zint2 = 0.0
 zint_t = 0.0
 zint_q = 0.0
 rhmean = 0.0
 zmean = 0.0
 do k=nlev,klev_max,-1
   if (icond(k) == 2) then
      cvgu = cvgu + max(0.0,ptenq(k))*(ph(k+1) - ph(k))/rg
      cvgh = cvgh + max(0.0,(ptent(k) + alpha_env(k)*ptenq(k)))*(ph(k+1) - ph(k))/rg
      zint = zint + (qvc(k) - qv(k) + (tc(k) - t(k))/alpha_env(k))*(ph(k+1) - ph(k))/rg
      zint2 = zint2 + (tc(k) - t(k) + alpha_env(k)*(qvc(k) - qv(k)))*(ph(k+1) - ph(k))/rg
      zint_t = zint_t + (tc(k) - t(k))*(ph(k+1) - ph(k))/rg
      zint_q = zint_q + (qvc(k) - qv(k))*(ph(k+1) - ph(k))/rg
      rhmean = rhmean + rh(k)*(ph(k+1) - ph(k))/rg
      zmean =  zmean + (ph(k+1) - ph(k))/rg
   endif
 enddo
 if (zmean > 0.0) then
   rhmean = rhmean/zmean
 else
   rhmean = 0.0
   print *,"Warning : no convective layer has been found ! RH mean set to zero"
 endif    
!
! Empirical moistening parameter for Kuo-Anthes scheme
! 
 bkuo = (1 - rhmean/1. - 0.1)**1
!print *,'**** kuo moistening parameter *****',bkuo
! 
!  T and q tendencies for Kuo (1965) convection scheme
! 
 rain = 0.0
 dtdt = 0.0
 dqdt = 0.0
 do k=nlev,klev_max,-1
   if (icond(k) == 2 .and. ptenq(k) > 0.0) then
     dtdt(k) = cvgu/zint*(tc(k) - t(k))*zcv1
     rain = rain + dtdt(k)*(ph(k+1) - ph(k))/(rg*alpha_env(k))
     dqdt(k) = (-ptenq(k) + cvgu/zint*(qvc(k) - qv(k)))*zcv2
     print *,'convective tendencies KUO',p(k),k,dtdt(k),dqdt(k),ptenq(k)*zcv2
   endif 
 enddo
 do k=nlev,1,-1   
   write (200,*) p(k),dtdt(k),dqdt(k)
 enddo
 print *,'surface rainfall',rain,zcape ! mm/jour
! write (150,*) imax,zlat,zlon,rain,zcape
! 
!  T and q tendencies for Kuo-Geleyn (1985) convection scheme
! 
 rain = 0.0
 dtdt = 0.0
 dqdt = 0.0
 do k=nlev,klev_max,-1
   if (icond(k) == 2 .and. ptenq(k) > 0.0) then
     dtdt(k) = cvgu/zint*((tc(k) - t(k)) + alpha_env(k)*(qvc(k) - qwb(k)))*zcv1
     rain = rain + dtdt(k)*(ph(k+1) - ph(k))/(rg*alpha_env(k))
     dqdt(k) = (-ptenq(k)+ cvgu/zint*(qwb(k) - qv(k)))*zcv2
     print *,'convective tendencies KUO-GELEYN',p(k),k,dtdt(k),dqdt(k),ptent(k)*zcv1
   endif 
   Cps = Cp*(1.0 - qv(k)) + qv(k)*Cpv
   write (300,*) p(k)/100.,t(k),twb(k),tc(k),qv(k),qwb(k),qvc(k),gz(k), &
               & Cps*t(k) + Lh(t(k),lsub)*qv(k) + gz(k)
 enddo
 do k=nlev,1,-1   
   write (201,*) p(k),dtdt(k),dqdt(k)
 enddo
 print *,'surface rainfall',rain,zcape ! mm/jour
! write (151,*) imax,zlat,zlon,rain,zcape
! 
!  T and q tendencies for Kuo symmetric convection scheme
! 
 rain = 0.0
 dtdt = 0.0
 dqdt = 0.0
 print *,'zint2=',zint2,'cvgh=',cvgh
 do k=nlev,klev_max,-1
   if (icond(k) == 2 .and. (ptent(k) + alpha_env(k)*ptenq(k)) > 0.0 .and. ptent(k) < 0.0) then
     dtdt(k) = (-ptent(k) + cvgh/zint2*(tc(k) - t(k)))*zcv1 
     rain = rain + dtdt(k)*(ph(k+1) - ph(k))/(rg*alpha_env(k))
     dqdt(k) = (-ptenq(k)+ cvgh/zint2*(qvc(k) - qv(k)))*zcv2
     print *,'convective tendencies KUO SYMMETRIC',p(k),k,dtdt(k),dqdt(k),ptent(k)*zcv1
   endif 
 enddo
 do k=nlev,1,-1   
   write (202,*) p(k),dtdt(k),dqdt(k)
 enddo
 print *,'surface rainfall',rain,zcape,cvgh, cvgu ! mm/jour
! write (152,*) imax,zlat,zlon,rain,zcape
! 
!  T and q tendencies for Kuo-Anthes 1977 convection scheme
! 
 rain = 0.0
 dtdt = 0.0
 dqdt = 0.0
 do k=nlev,klev_max,-1
   if (icond(k) == 2 .and. ptenq(k) > 0.0) then
     dtdt(k) = alpha_env(k)*(1.- bkuo)*cvgu/zint_t*(tc(k) - t(k))*zcv1 
     rain = rain + dtdt(k)*(ph(k+1) - ph(k))/(rg*alpha_env(k))
     dqdt(k) = (-ptenq(k)+ bkuo*cvgu/zint_q*(qvc(k) - qv(k)))*zcv2
     print *,'convective tendencies KUO-ANTHES',p(k),k,dtdt(k),dqdt(k),ptent(k)*zcv1
   endif 
 enddo
 do k=nlev,1,-1   
   write (203,*) p(k),dtdt(k),dqdt(k)
 enddo
 print *,'surface rainfall',rain,zcape,cvgh, cvgu ! mm/jour
! write (153,*) imax,zlat,zlon,rain,zcape 
!
 close (unit=200)
 close (unit=201)
 close (unit=202)
 close (unit=203)
 close (unit=300)
!
 stop
end program kuo_schemes
