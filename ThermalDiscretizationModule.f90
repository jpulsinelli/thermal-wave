MODULE ThermalDiscretizationModule

  USE ThermalParametersModule, ONLY: &
       DP
  USE OMP_LIB
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Heat

CONTAINS


  SUBROUTINE ComputeIncrement_Heat &
       ( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, X_L, X_R, dX, dY, dZ, dt, time, del, T, dU )

    INTEGER,  INTENT(in)     :: nE, nX, nY, nZ
    INTEGER,  INTENT(in)     :: iE_L, iX_L, iY_L, iZ_L
    INTEGER,  INTENT(in)     :: iE_R, iX_R, iY_R, iZ_R
    REAL(DP), INTENT(in)    :: X_L, X_R, del, time
    REAL(DP), INTENT(in)    :: dX, dY, dZ, dt
    REAL(DP), INTENT(inout) :: T (1:nE,0:nX+1,0:nY+1,0:nZ+1)
    REAL(DP), INTENT(out)   :: dU(1:nE,1:nX+0,1:nY+0,1:nZ+0)

    INTEGER :: iE, iX, iY, iZ
    REAL(DP) :: T_star(1:nE,1:nX,1:nY,1:nZ)
    REAL(DP) :: D_T(1:nE,1:nX,1:nY,1:nZ), U(1:nE,1:nX,1:nY,1:nZ)
    REAL(DP) :: F(1:nE,1:nX,1:nY,1:nZ), F_prime(1:nE,1:nX,1:nY,1:nZ)
    REAL(DP) :: analytic(1:nE,1:nX,1:nY,1:nZ), diff(1:nE,1:nX,1:nY,1:nZ)
    REAL(DP) :: F_zero, rtol, atol

      rtol = 1.d-7
      atol = 1.d-8
    
      CALL ApplyBoundaryConditions( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, X_L, X_R, dX, T, time )

! $OMP TARGET DATA MAP(to:T) MAP(from:dU)  !GPU only
! $OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) !GPU only
! $OMP PARALLEL DO COLLAPSE(4) !CPU only
    DO iZ = iZ_L, iZ_R
      DO iY = iY_L, iY_R
        DO iX = iX_L, iX_R
          DO iE = iE_L, iE_R
              
             D_T(iE,iX,iY,iZ) = (1/(dx*dx))*(T(iE,iX-1,iY,iZ) - 2*T(iE,iX,iY,iZ) + T(iE,iX+1,iY,iZ)) &
                + (1/(dy*dy))*(T(iE,iX,iY-1,iZ) - 2*T(iE,iX,iY,iZ) + T(iE,iX,iY+1,iZ)) &
                + (1/(dz*dz))*(T(iE,iX,iY,iZ-1) - 2*T(iE,iX,iY,iZ) + T(iE,iX,iY,iZ+1))
             T_star(iE,iX,iY,iZ) = T(iE,iX,iY,iZ) + dt*D_T(iE,iX,iY,iZ)
             U(iE,iX,iY,iZ) = T_star(iE,iX,iY,iZ)
             F(iE,iX,iY,iZ) = U(iE,iX,iY,iZ) - 8*(dt/del)*U(iE,iX,iY,iZ)**2 + 8*(dt/del)*U(iE,iX,iY,iZ)**3 - T_star(iE,iX,iY,iZ)
             F_zero = abs(F(iE,iX,iY,iZ))

             DO WHILE (abs(F(iE,iX,iY,iZ)) > rtol*F_zero + atol)
                F_prime(iE,iX,iY,iZ) = 1 - 16*(dt/del)*U(iE,iX,iY,iZ) + 24*(dt/del)*U(iE,iX,iY,iZ)**2
                diff(iE,iX,iY,iZ) = -F(iE,iX,iY,iZ)/F_prime(iE,iX,iY,iZ)
                U(iE,iX,iY,iZ) = U(iE,iX,iY,iZ) + diff(iE,iX,iY,iZ)
                F(iE,iX,iY,iZ) = U(iE,iX,iY,iZ) - 8*(dt/del)*U(iE,iX,iY,iZ)**2 + 8*(dt/del)*U(iE,iX,iY,iZ)**3 - T_star(iE,iX,iY,iZ)
             ENDDO

             dU(iE,iX,iY,iZ) = U(iE,iX,iY,iZ) - T(iE,iX,iY,iZ)
          
          ENDDO
        ENDDO
      ENDDO
    ENDDO
! $OMP END PARALLEL DO
! $OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
! $OMP END TARGET DATA
    
  END SUBROUTINE ComputeIncrement_Heat


  SUBROUTINE ApplyBoundaryConditions &
    ( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, X_L, X_R, dX, U, time )

    INTEGER,  INTENT(in)    :: nE, nX, nY, nZ
    INTEGER,  INTENT(in)    :: iE_L, iX_L, iY_L, iZ_L
    INTEGER,  INTENT(in)    :: iE_R, iX_R, iY_R, iZ_R
    REAL(DP), INTENT(in)    :: X_L, X_R, dX, time
    REAL(DP), INTENT(inout) :: U(1:nE,0:nX+1,0:nY+1,0:nZ+1)

    ! --- X-Dimension
    
      IF (iX_L == 1) THEN
         U(iE_L:iE_R, 0, iY_L:iY_R, iZ_L:iZ_R) =  0.5 - 0.5*tanh(X_L - 0.5*dX- 2*time)
      ENDIF

      IF (iX_R == nX) THEN
         U(iE_L:iE_R, nX+1, iY_L:iY_R, iZ_L:iZ_R) = 0.5 - 0.5*tanh(X_R + 0.5*dX- 2*time)
      ENDIF

    ! --- Y-Dimension
      
       IF (iY_L == 1) THEN
         U(iE_L:iE_R, iX_L:iX_R, 0, iZ_L:iZ_R) = U(iE_L:iE_R, iX_L:iX_R, nY, iZ_L:iZ_R)
      ENDIF

      IF (iY_R == nY) THEN
         U(iE_L:iE_R, iX_L:iX_R, nY+1, iZ_L:iZ_R) = U(iE_L:iE_R, iX_L:iX_R, 1, iZ_L:iZ_R)
      ENDIF

    ! --- Z-Dimension

      IF (iZ_L == 1) THEN
         U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, 0) = U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, nZ)
      ENDIF

      IF (iZ_R == nZ) THEN
         U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, nZ+1) = U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, 1)
      ENDIF
    
    
    
  END SUBROUTINE ApplyBoundaryConditions

  
END MODULE ThermalDiscretizationModule
