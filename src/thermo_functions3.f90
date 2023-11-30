!----------------------------------------
 FUNCTION esat(t)
!************************************************************
! saturation water vapour pressure
!
! input : t in K
! ouput : esat in Pa
!
! Jean-Francois Mahfouf (30/11/2023)
!
! Formula proposed by Huang (2018) J.A.M.C., 57, 1265-1272
!
!************************************************************
 IMPLICIT NONE
 REAL :: esat
 REAL, INTENT(IN) :: t
 REAL, PARAMETER :: aw=34.494, bw=4924.99, cw=1.57, d1w=237.1, d2w=105.0
 REAL, PARAMETER :: ai=43.494, bi=6546.80, ci=2.00, d1i=278.0, d2i=868.0 
 REAL :: t_celsius
 t_celsius = t - 273.15
 IF (t >= 273.15) THEN
   esat=exp(aw - bw/(t_celsius + d1w))/(t_celsius + d2w)**cw 
 ELSE
   esat=exp(ai - bi/(t_celsius + d1i))/(t_celsius + d2i)**ci
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
! Jean-Francois Mahfouf (30/11/2023)
!
! Formula proposed by Huang (2018) J.A.M.C., 57, 1265-1272
!
!************************************************************
 IMPLICIT NONE
 REAL :: desat, esat
 REAL, INTENT(IN) :: t
 REAL, PARAMETER :: aw=34.494, bw=4924.99, cw=1.57, d1w=237.1, d2w=105.0
 REAL, PARAMETER :: ai=43.494, bi=6546.80, ci=2.00, d1i=278.0, d2i=868.0 
 REAL :: t_celsius 
 t_celsius = t - 273.15
 IF (t >= 273.15) THEN    
   esat=exp(aw - bw/(t_celsius + d1w))/(t_celsius + d2w)**cw      
   desat=esat*(bw/(t_celsius + d1w)**2 - cw/(t_celsius + d2w))
 ELSE
   esat=exp(ai - bi/(t_celsius + d1i))/(t_celsius + d2i)**ci    
   desat=esat*(bi/(t_celsius + d1i)**2 - ci/(t_celsius + d2i))      
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
