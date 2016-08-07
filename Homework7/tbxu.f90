! complie: gfortran -o tbxu tbxu.f90 -lblas -llapack

PROGRAM tbxu
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER :: rcut = 2.6D0   ! Cutoff distance
   DOUBLE PRECISION, PARAMETER :: infinity = HUGE(1.0D0)
   DOUBLE PRECISION , ALLOCATABLE :: bondcosines(:,:,:)
   DOUBLE PRECISION , ALLOCATABLE :: bondlength(:,:)
   DOUBLE PRECISION :: ebs
   DOUBLE PRECISION :: erep
   DOUBLE PRECISION :: etot
   DOUBLE PRECISION :: etotperatom
   DOUBLE PRECISION :: alat(1:3,1:3)
   INTEGER :: natoms
   DOUBLE PRECISION , ALLOCATABLE :: tau(:,:)
   WRITE(*,*)
   WRITE(*,"(A)") "Starting program TBXU..."
   WRITE(*,*)
   CALL read_alat(alat)
   CALL read_number_of_atoms(natoms)
   ALLOCATE(tau(1:natoms,1:3))
   CALL read_tau(tau)
   ALLOCATE(bondlength(1:natoms,1:natoms))
   ALLOCATE(bondcosines(1:natoms,1:natoms,1:3))
   CALL compute_bond_info(alat, tau, bondlength, bondcosines)
   CALL evaluate_erep(bondlength, erep)
   CALL evaluate_ebs(bondlength, bondcosines, ebs)
   etot = ebs + erep
   WRITE(*,"(A)") "Adding up energy contributions..."
   WRITE(*,"(3X,A,F16.10)") "TOTAL ENERGY:          ", etot
   etotperatom = etot/dble(natoms)
   WRITE(*,"(3X,A,F16.10)") "TOTAL ENERGY PER ATOM: ", etotperatom
   WRITE(*,*)
   DEALLOCATE(bondcosines)
   DEALLOCATE(bondlength)
   DEALLOCATE(tau)
   WRITE(*,"(A)") "Finishing program TBXU..."
   WRITE(*,*)
CONTAINS

   !----------------------------------------------------------------------------
   SUBROUTINE read_alat(alat)
      ! This subroutine reads from standard input the lattice vectors of a
      ! crystallographic structure, and returns them in a 3x3 array alat.
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: alat(:,:)
                            INTEGER :: i
   WRITE(*,"(A)") "Reading lattice vectors..."
   DO i = 1, 3
      READ(*,*) alat(i,1:3)
      WRITE(*,"(3X,3F16.10)") alat(i,1:3)
   END DO
   WRITE(*,*)
END SUBROUTINE read_alat


!----------------------------------------------------------------------------
SUBROUTINE read_number_of_atoms(n)
   ! This subroutine reads from standard input the number of atoms to
   ! be considered in the simulation cell, and returns it.
   IMPLICIT NONE
               INTEGER, INTENT(OUT) :: n
   WRITE(*,"(A)") "Reading number of atoms..."
   READ(*,*) n
   WRITE(*,"(3X,I6)") n
   WRITE(*,*)
END SUBROUTINE read_number_of_atoms
!----------------------------------------------------------------------------


SUBROUTINE read_tau(tau)
   ! This subroutine reads from standard input the reduced coordinates
   ! of the atoms in the cell, and returns them in an array.
   IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: tau(:,:)
                            INTEGER :: i
                            INTEGER :: n
   WRITE(*,"(A)") "Reading reduced coordinates..."
   n = size(tau,1)
   DO i = 1, n
      READ(*,*) tau(i,1:3)
      WRITE(*,"(3X,3F16.10)") tau(i,1:3)
   END DO
   WRITE(*,*)
   
END SUBROUTINE read_tau

!----------------------------------------------------------------------------
SUBROUTINE compute_bond_info(alat, tau, length, cosines)
   ! Given the lattice vectors of a cell (alat) and the reduced positions
   ! of a set of atoms in this cell (tau), this subroutine computes the
   ! distances (length) and direction cosines (cosines) for every
   ! pair of atoms, taking into account periodic boundary conditions
   ! (i.e., the distance is the shortest between one atom and any of the
   ! periodic copies of the other).
   IMPLICIT NONE
       DOUBLE PRECISION, INTENT(IN) :: alat(:,:)
       DOUBLE PRECISION, INTENT(IN) :: tau(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: cosines(:,:,:)
      DOUBLE PRECISION, INTENT(OUT) :: length(:,:)
                   DOUBLE PRECISION :: distsq
                            INTEGER :: i
                            INTEGER :: j
                            INTEGER :: kcell1
                            INTEGER :: kcell2
                            INTEGER :: kcell3
                   DOUBLE PRECISION :: mindistsq
                            INTEGER :: n
                   DOUBLE PRECISION :: newposj(1:3)
     DOUBLE PRECISION , ALLOCATABLE :: pos(:,:)
                   DOUBLE PRECISION :: rcutsq
                   DOUBLE PRECISION :: shortestbond
   WRITE(*,"(A)") "Computing atomic distances and direction cosines..."
   n = SIZE(tau,1)
   ALLOCATE(pos(1:n,1:3))
   DO i = 1, n
      pos(i,1) = tau(i,1)*alat(1,1) + tau(i,1)*alat(2,1) + tau(i,1)*alat(3,1)
      pos(i,2) = tau(i,2)*alat(1,2) + tau(i,2)*alat(2,2) + tau(i,2)*alat(3,2)
      pos(i,3) = tau(i,3)*alat(1,3) + tau(i,3)*alat(2,3) + tau(i,3)*alat(3,3)
END DO
   shortestbond = rcut
   rcutsq = rcut*rcut
   DO i=1, n-1
      DO j=i+1,n
         mindistsq = infinity
         DO kcell1=-1,1
            DO kcell2=-1,1
               DO kcell3=-1,1
                  newposj(1) = pos(j,1) + kcell1*alat(1,1)  &
                                        + kcell2*alat(2,1)  &
                                        + kcell3*alat(3,1)
                  newposj(2) = pos(j,2) + kcell1*alat(1,2)  &
                                        + kcell2*alat(2,2)  &
                                        + kcell3*alat(3,2)
                  newposj(3) = pos(j,3) + kcell1*alat(1,3)  &
                                        + kcell2*alat(2,3)  &
                                        + kcell3*alat(3,3)
                  distsq = ( pos(i,1) - newposj(1) )**2 +  &
                           ( pos(i,2) - newposj(2) )**2 +  &
                           ( pos(i,3) - newposj(3) )**2
                  IF (distsq .LE. mindistsq) THEN
                     mindistsq = distsq
                     length(i,j) = sqrt(distsq)
                     cosines(i,j,1) = (newposj(1)-pos(i,1))/sqrt(distsq)
                     cosines(i,j,2) = (newposj(2)-pos(i,2))/sqrt(distsq)
                     cosines(i,j,3) = (newposj(3)-pos(i,3))/sqrt(distsq)
                     length(j,i) = length(i,j)
                     cosines(j,i,1) = -cosines(i,j,1)
                     cosines(j,i,2) = -cosines(i,j,2)
                     cosines(j,i,3) = -cosines(i,j,3)

                  END IF 
              END DO
          END DO 
        END DO
         IF (length(i,j) < shortestbond) THEN
            shortestbond = length(i,j)
         END IF 
    END DO
  END DO
   DEALLOCATE(pos)
   WRITE(*,"(3X,A,F16.10)") "Shortest bond:", shortestbond
   WRITE(*,*)
END SUBROUTINE compute_bond_info

!----------------------------------------------------------------------------
SUBROUTINE evaluate_erep(bondlength, erep)
   ! Given a collection of atoms where for every pair we know
   ! the distance between them (bondlength(i,j) for atoms i and j), this
   ! subroutine computes the repulsion energy according to the paper
   ! "A transferable tight-binding potential for carbon", by Xu et al.
   IMPLICIT NONE
       DOUBLE PRECISION, INTENT(IN) :: bondlength(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: erep
        DOUBLE PRECISION , PARAMETER :: dm =  rcut

        DOUBLE PRECISION , PARAMETER :: phi0 =  8.18555D0
        DOUBLE PRECISION , PARAMETER :: m =  3.30304D0
        DOUBLE PRECISION , PARAMETER :: mc =  8.6655D0
        DOUBLE PRECISION , PARAMETER :: dc =  2.1052D0
        DOUBLE PRECISION , PARAMETER :: d0 =  1.64d0
        DOUBLE PRECISION , PARAMETER :: d1 =  2.57D0
        DOUBLE PRECISION , PARAMETER :: c0tphi =  2.2504290109D-8
        DOUBLE PRECISION , PARAMETER :: c1tphi = -1.4408640561D-6
        DOUBLE PRECISION , PARAMETER :: c2tphi =  2.1043303374D-5
        DOUBLE PRECISION , PARAMETER :: c3tphi =  6.6024390226D-5
        DOUBLE PRECISION , PARAMETER :: c0f    = -2.5909765118191D0

        DOUBLE PRECISION , PARAMETER :: c1f=  0.5721151498619D0

        DOUBLE PRECISION , PARAMETER :: c2f= -1.7896349903996D-3

        DOUBLE PRECISION , PARAMETER :: c3f=  2.3539221516757D-5

        DOUBLE PRECISION , PARAMETER :: c4f= -1.24251169551587D-7

                DOUBLE PRECISION :: f
                DOUBLE PRECISION :: erepperatom
                         INTEGER :: i
                         INTEGER :: j
                         INTEGER :: n
                DOUBLE PRECISION :: phi
                DOUBLE PRECISION :: r
                DOUBLE PRECISION :: rdiff
                DOUBLE PRECISION :: sumphi
                
                WRITE(*,"(A)") "Evaluating repulsion energy..."

                n = SIZE(bondlength,1)
                erep = 0.0D0
                DO i=1, n
                   sumphi = 0.0D0
                   DO j=1, n
                      IF (i .NE. j) THEN
                         r = abs(bondlength(i,j))
                         phi = 0.0
                IF (r>0.0D0) THEN
                    IF (r<d1) THEN
                        phi = phi0 * (d0/r)**m * EXP(m*(-(r/dc)**mc+(d0/dc)**mc))
                    ELSE IF ( r < dm ) THEN
                      rdiff = r - d1
                      phi = c0tphi + c1tphi*rdiff + c2tphi*rdiff*rdiff  &
                            + c3tphi*rdiff*rdiff*rdiff
                    END IF
                    sumphi = sumphi + phi
                END IF
            END IF
        END DO
        f = c0f + c1f*sumphi + c2f*sumphi*sumphi + c3f*sumphi*sumphi*sumphi  &
              + c4f*sumphi*sumphi*sumphi*sumphi
      erep = erep + f
   END DO
   WRITE(*,"(3X,A,F16.10)") "REPULSION ENERGY:          ", erep
   erepperatom = erep/n ;
   WRITE(*,"(3X,A,F16.10)") "REPULSION ENERGY PER ATOM: ", erepperatom
   WRITE(*,*)
END SUBROUTINE evaluate_erep
!----------------------------------------------------------------------------
SUBROUTINE evaluate_ebs(length, cosines, ebs)
   ! Given a collection of atoms in a line where for every pair we
   ! know the distance between them (bondlength(i,j) for atoms i and j),
   ! this subroutine computes the bandstructure energy according to the
   ! paper "A transferable tight-binding potential for carbon", by Xu et al.

IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: length(:,:)
    DOUBLE PRECISION, INTENT(IN) :: cosines(:,:,:)
   DOUBLE PRECISION, INTENT(OUT) :: ebs
             INTEGER , PARAMETER :: norb = 4
     DOUBLE PRECISION, PARAMETER :: es       = -2.99D0
     DOUBLE PRECISION, PARAMETER :: ep       =  3.71D0
     DOUBLE PRECISION, PARAMETER :: vsssigma = -5.0D0
     DOUBLE PRECISION, PARAMETER :: vspsigma =  4.7D0
     DOUBLE PRECISION, PARAMETER :: vppsigma =  5.5D0
    DOUBLE PRECISION, PARAMETER :: vpppi = -1.55D0
    DOUBLE PRECISION, PARAMETER :: rm=  rcut
    DOUBLE PRECISION, PARAMETER :: n =  2.0D0
    DOUBLE PRECISION, PARAMETER :: nc =  6.5D0
    DOUBLE PRECISION, PARAMETER :: rc =  2.18D0
    DOUBLE PRECISION, PARAMETER :: r0=  1.536329D0
    DOUBLE PRECISION, PARAMETER :: r1 =  2.45D0
    DOUBLE PRECISION, PARAMETER :: c0ts=  6.7392620074314D-3
    DOUBLE PRECISION, PARAMETER :: c1ts= -8.1885359517898D-2
    DOUBLE PRECISION, PARAMETER :: c2ts=  0.1932365259144D0
    DOUBLE PRECISION, PARAMETER :: c3ts=  0.3542874332380D0

    DOUBLE PRECISION :: ebsperatom
  DOUBLE PRECISION , ALLOCATABLE :: eigval(:)
    DOUBLE PRECISION :: ess
    DOUBLE PRECISION :: esx
    DOUBLE PRECISION :: esy
    DOUBLE PRECISION :: esz
    DOUBLE PRECISION :: exx
    DOUBLE PRECISION :: exy
    DOUBLE PRECISION :: exz
    DOUBLE PRECISION :: eyy
    DOUBLE PRECISION :: eyz
    DOUBLE PRECISION :: ezz
  DOUBLE PRECISION , ALLOCATABLE :: hamiltonian(:,:)
     INTEGER :: i
     INTEGER :: ih
     INTEGER :: info
     INTEGER :: j
     INTEGER :: jh
    DOUBLE PRECISION :: ll
     INTEGER :: lwork
    DOUBLE PRECISION :: mm
     INTEGER :: natoms
     INTEGER :: nh
    DOUBLE PRECISION :: nn
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: rdiff
    DOUBLE PRECISION :: sr
    DOUBLE PRECISION , ALLOCATABLE :: work(:)
        WRITE(*,"(A)") "Evaluating band-structure energy..."

      natoms = SIZE(length,1)
      nh = natoms*norb
      ALLOCATE(hamiltonian(1:nh,1:nh))
      hamiltonian(:,:) = 0.0D0
!     Compute value of diagonal terms
      DO i = 1, natoms
         ih = norb*(i-1)
         hamiltonian(ih+1,ih+1) = es
         hamiltonian(ih+2,ih+2) = ep
         hamiltonian(ih+3,ih+3) = ep
         hamiltonian(ih+4,ih+4) = ep
    END DO
    !     Compute value of off-diagonal terms (upper part of the matrix)
          DO i = 1, natoms-1
             ih = norb*(i-1)
             DO j = i+1, natoms
                jh = norb*(j-1)
                r = abs(length(i,j)) 
                sr = 0.0D0
                IF (r>0.0)THEN
                    IF ( r < r1 ) THEN
                        sr = (r0/r)**n * EXP(n*(-(r/rc)**nc+(r0/rc)**nc))
                    ELSE IF ( r < rm ) THEN
                        rdiff = r - r1
                        sr = c0ts + c1ts*rdiff + c2ts*rdiff*rdiff  &
                            + c3ts*rdiff*rdiff*rdiff
                    END IF
               ll = cosines(i,j,1)
               mm = cosines(i,j,2)
               nn = cosines(i,j,3)
               ess = vsssigma
               esx = ll * vspsigma
               esy = mm * vspsigma
               esz = nn * vspsigma
               exx = ll**2 * vppsigma + (1-ll**2) * vpppi
               eyy = mm**2 * vppsigma + (1-mm**2) * vpppi
               ezz = nn**2 * vppsigma + (1-nn**2) * vpppi
               exy = ll*mm * (vppsigma-vpppi)
               eyz = mm*nn * (vppsigma-vpppi)
               exz = ll*nn * (vppsigma-vpppi)
               hamiltonian(ih+1,jh+1) =  ess * sr
               hamiltonian(ih+2,jh+2) =  exx * sr
               hamiltonian(ih+3,jh+3) =  eyy * sr
               hamiltonian(ih+4,jh+4) =  ezz * sr
               hamiltonian(ih+1,jh+2) = -esx * sr
               hamiltonian(ih+1,jh+3) = -esy * sr
               hamiltonian(ih+1,jh+4) = -esz * sr
               hamiltonian(ih+2,jh+1) =  esx * sr
               hamiltonian(ih+3,jh+1) =  esy * sr
               hamiltonian(ih+4,jh+1) =  esz * sr
               hamiltonian(ih+2,jh+3) =  exy * sr
               hamiltonian(ih+2,jh+4) =  exz * sr
               hamiltonian(ih+3,jh+4) =  eyz * sr
               hamiltonian(ih+3,jh+2) =  exy * sr
               hamiltonian(ih+4,jh+2) =  exz * sr
               hamiltonian(ih+4,jh+3) =  eyz * sr
            END IF
        END DO 
    END DO
!     Diagonalize Hamiltonian (the upper part of the matrix is given):
      lwork = MAX(1,3*nh-1)
      ALLOCATE(work(1:lwork))
      ALLOCATE(eigval(1:nh))
!      CALL dsyev("N","U",nh,hamiltonian,nh,eigval,work,lwork,info)
      CALL dsyev('N','U',nh,hamiltonian,nh,eigval,work,lwork,info)
      WRITE(*,"(3X,A,I3)") "Diagonalization performed, info equals ",info
      WRITE(*,"(3X,A)") "EIGENVALUES:"
      WRITE(*,"(4F18.6)") eigval(:)
      ebs = 0.0
      DO i=1, nh/2
         ebs = ebs + eigval(i)
      END DO
      ebs = ebs*2
      
      WRITE(*,"(3X,A,F16.10)") "BAND-STRUCTURE ENERGY:"
      ebsperatom = ebs/dble(natoms) ;
      WRITE(*,"(3X,A,F16.10)") "BAND-STRUCTURE ENERGY PER ATOM: ", ebsperatom
      WRITE(*,*)
   END SUBROUTINE evaluate_ebs
   !----------------------------------------------------------------------------
END PROGRAM