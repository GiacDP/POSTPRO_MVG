 SUBROUTINE welch(fsig,fspec,m,n,kseg)
 USE fftw3
!
 INTEGER, INTENT(IN) :: m,n,kseg
 REAL, DIMENSION(n), INTENT(IN) :: fsig
 REAL, DIMENSION(m+1), INTENT(OUT) :: fspec
 TYPE(C_PTR) :: plan
 REAL(C_DOUBLE), dimension(2*m) :: datare,wind
 COMPLEX(C_DOUBLE_COMPLEX), dimension(m+1) :: datacmplx
!
 REAL :: wss
!
 pi = atan(1.)*4.
 m2 = m+m
 fspec = 0. 
 loff  = 0
 fm  = float(m)
 fm2 = float(m2)
 wss = 0.
!
 DO l=1,m2
  wind(l) = 1.                              ! Rectangular
  ! wind(l) = 0.5*(1.-cos(2.*pi*(l-1)/fm2)) ! Hann
  ! wind(l) = 1.-(((float(l)-1.)-fm)/fm)**2 ! Welch
  wss = wss+wind(l)*wind(l)
 ENDDO
!
 plan  = fftw_plan_dft_r2c_1d(m2,datare,datacmplx,FFTW_ESTIMATE)
!
 DO k=1,kseg
  datare = fsig(1+loff:m2+loff)
  ! Remove mean
  ! datarem = 0.
  ! DO l=1,m2
  !  datarem = datarem+datare(l)
  ! ENDDO
  ! datare = datare-datarem/fm2
  DO l=1,m2
   datare(l) = datare(l)*wind(l)
  ENDDO
  CALL fftw_execute_dft_r2c(plan,datare,datacmplx)
  fspec  = fspec+(real(datacmplx)**2+aimag(datacmplx)**2)
  loff   = loff+m
 ENDDO
!
 fspec = fspec/(kseg*wss)
!
 CALL fftw_destroy_plan(plan)
!
 RETURN
 END SUBROUTINE welch
