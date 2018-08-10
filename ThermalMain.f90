  PROGRAM ThermalMain

  USE ThermalParametersModule, ONLY: &
       DP, nE, nX, nY, nZ, Pi, &
       X_L, X_R, Y_L, Y_R, Z_L, Z_R, t_end, C, del
  USE ThermalDiscretizationModule, ONLY: &
       ComputeIncrement_Heat

  USE OMP_LIB
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: iE, iX, iY, iZ, mpierr, iCycle, threads
  REAL(DP) :: dX, dY, dZ, t, dt, wTime, err, t_cell
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: U, dU, analytic, diff, Ureal

  CALL MPI_INIT( mpierr )
  PRINT*, "mpierr = ", mpierr

  dX = ( X_R - X_L ) / DBLE( nX )
  dY = ( Y_R - Y_L ) / DBLE( nY )
  dZ = ( Z_R - Z_L ) / DBLE( nZ )
  
  ALLOCATE( U (1:nE,0:nX+1,0:nY+1,0:nZ+1) )
  ALLOCATE( dU(1:nE,1:nX+0,1:nY+0,1:nZ+0) )
  ALLOCATE( analytic(1:nE,1:nX,1:nY,1:nZ) )
  ALLOCATE( Ureal(1:nE,1:nX,1:nY,1:nZ) )
  ALLOCATE( diff(1:nE,1:nX,1:nY,1:nZ) )
  
  wTime = MPI_WTIME( )


CALL OMP_SET_NUM_THREADS(16) !Unnecessary for GPU, helps with CPU
!The following omp lines may mildly reduce CPU performance
! $OMP TARGET DATA MAP(from:U)
! $OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4)
  DO iZ = 1, nZ
    DO iY = 1, nY
      DO iX = 1, nX
        DO iE = 1, nE
            
           U(iE,iX,iY,iZ) = 0.5 - 0.5*tanh(X_L + (DBLE(iX)-0.5)*dX)
           threads = OMP_GET_NUM_THREADS()
            
        ENDDO
      ENDDO
    ENDDO
 ENDDO
! $OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
! $OMP END TARGET DATA 

  wTime = MPI_WTIME( ) - wTime
write(*,*) threads

  PRINT*, "wTime (allocation loop) = ", wTime
 
  t = 0.0_DP
  dt= C*dx*dx
 
  wTime = MPI_WTIME( )

  iCycle = 0

  DO WHILE( t < t_end )
  !DO WHILE(iCycle < 100)

    iCycle = iCycle + 1

    CALL ComputeIncrement_Heat &
         ( nE, nX, nY, nZ, 1, nE, 1, nX, 1, nY, 1, nZ, X_L, X_R, dX, dY, dZ, dt, t, del, U, dU )

    U(1:nE,1:nX,1:nY,1:nZ) = U(1:nE,1:nX,1:nY,1:nZ) + dU(1:nE,1:nX,1:nY,1:nZ)
    
    t = t + dt

  END DO


  wTime = MPI_WTIME( ) - wTime
 
  PRINT*, "wTime (computational loop) = ", wTime

 DO iZ = 1, nZ
    DO iY = 1, nY
      DO iX = 1, nX
        DO iE = 1, nE
           analytic(iE,iX,iY,iZ) = 0.5*(1-tanh(X_L + (DBLE(iX)-0.5)*dX - 2*t))
        ENDDO
      ENDDO
    ENDDO
 ENDDO

  Ureal = U(1:nE,1:nX,1:nY,1:nz)
  diff = abs(Ureal - analytic)
  err = maxval(diff)
 
  OPEN(unit=10, file = 'ThermalFinal.dat')
  WRITE(10,*) U(1,:,1,1)
  CLOSE(10)
  DEALLOCATE( U, dU, analytic, diff, Ureal )
  t_cell=wTime/nX/nY/nZ/nE/iCycle
  
  CALL MPI_FINALIZE( mpierr )

  WRITE(*,*) 'Advanced to time',t,'in',iCycle, 'steps'
  WRITE(*,*) 'Time per cell', t_cell
  WRITE(*,*) 'nE=',nE
  WRITE(*,*) 'nX=',nX
  WRITE(*,*) 'nY=',nY 
  WRITE(*,*) 'nZ=',nZ
  WRITE(*,*) 'dt=',dt
  WRITE(*,*) 'error=',err

END PROGRAM ThermalMain

