!----------------------------------------
 FUNCTION esat(t)
!****************************************
! saturation water vapour pressure
!
! input : t in K
! ouput : esat in Pa
!
! Jean-Francois Mahfouf (03/07/2001)
!
!****************************************
 IMPLICIT NONE
 REAL :: esat
 REAL, INTENT(IN) :: t
 REAL, PARAMETER :: aw=611.21, bw=17.502, cw=32.19, t00=273.16
 REAL, PARAMETER :: ai=611.21, bi=22.587, ci=-0.7
 IF (t >= t00) THEN
   esat=aw*exp(bw*(t - t00)/(t - cw)) ! ECMWF
 ELSE
   esat=ai*exp(bi*(t - t00)/(t - ci)) ! ECMWF
 ENDIF
 END FUNCTION esat
!---------------------------------------
 FUNCTION qsat(p,t)
!****************************************
! specific humidity at saturation
!
! inputs : p in Pa
!          t in K
!  ouput : qsat in kg/kg
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*****************************************
 IMPLICIT NONE
 REAL :: qsat,eps,esat
 REAL, INTENT(IN) :: p,t
 REAL, PARAMETER :: rv=461.5,rd=287.04
 eps=rd/rv
 qsat=eps*esat(t)/(p-esat(t)*(1.0-eps))
 END FUNCTION qsat
!------------------------------------------------------------
 FUNCTION desat(t)
!************************************************************
! derivative of saturation water vapour pressure w.r.t.to t
!
! input : t in K
! ouput : desat in Pa/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!************************************************************
 IMPLICIT NONE
 REAL :: desat
 REAL, INTENT(IN) :: t
 REAL, PARAMETER :: aw=611.21, bw=17.502, cw=32.19, t00=273.16
 REAL, PARAMETER :: ai=611.21, bi=22.587, ci=-0.7
 IF (t >= t00) THEN
   desat=aw*bw*(t00 - cw)*exp(bw*(t - t00)/(t - cw))/ &
&        ((t - cw)*(t - cw))          
 ELSE
   desat=ai*bi*(t00 - ci)*exp(bi*(t-t00)/(t - ci))/ &
&        ((t - ci)*(t - ci))          
 ENDIF
 END FUNCTION desat 
!--------------------------------------------------------
 FUNCTION dqsat(p,t)
!********************************************************
! derivative of specific humidity at saturation w.r.t t
!
! inputs : p in Pa
!          t in K
!  ouput : qsat in kg/kg/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*********************************************************
 IMPLICIT NONE
 REAL :: dqsat,eps,desat,esat
 REAL, INTENT(IN) :: p,t
 REAL, PARAMETER :: rv=461.5,rd=287.04
 eps=rd/rv
 dqsat=eps*desat(t)*p/(p-esat(t)*(1.0-eps))**2
 END FUNCTION dqsat
!--------------------------------------------------------
 FUNCTION Lh(t,lsub)
!********************************************************
! Latent heat of vaporisation or sublimation 
!
! inputs : t in K
!          lsub is a logical for sublimation
!
! ouput : Lh in J/(kg/kg)
!
! Jean-Francois Mahfouf (24/11/2023)
!
!*********************************************************
 IMPLICIT NONE
 REAL :: Lh, Lv, Li
 REAL, INTENT(IN) :: t
 LOGICAL, INTENT(IN) :: lsub
 REAL, PARAMETER :: Lv0=2.501E6, Li0=0.3337E6
 REAL, PARAMETER :: Cpv=1846.0, Cl=4218.0, Ci=2016.0, t00=273.16
 Lv = Lv0 + (t - t00)*(Cpv - Cl)
 Li = Li0 + (t - t00)*(Cl - Ci)
 Lh = Lv
 if (t < t00 .and. lsub) Lh = (Li + Lv) 
 END FUNCTION Lh
