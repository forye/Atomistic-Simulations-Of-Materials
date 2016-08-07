! gfortran -o  eigen_vals eigenvalues.f90 -llapack

PROGRAM eigenvalues
    DOUBLE PRECISION, ALLOCATABLE :: eigval(:)
    INTEGER :: i
    INTEGER :: info
    INTEGER :: lwork
    DOUBLE PRECISION, ALLOCATABLE :: matrix(:,:)
    INTEGER :: n
    DOUBLE PRECISION, ALLOCATABLE :: work(:)
    
    ! Create matrix for this example:
    n = 32
    ALLOCATE(matrix(1:n,1:n))
    DO i = 1, n
       matrix(i,i) = 2.0d0
    END DO
    
    DO i = 1, n-1
       matrix(i,i+1) = -1.0d0
    END DO
    DO i = 2, n
       matrix(i-1,i) = -1.0d0
    END DO
    
    ! Call library to diagonalize matrix:
    lwork = MAX(1,3*n-1)    
    ALLOCATE(work(1:lwork))
    ALLOCATE(eigval(1:n))
    
    CALL dsyev('N','U',n,matrix,n,eigval,work,lwork,info)
    
    WRITE(*,"(3X,A,I3)") 'Diagonalization performed, info equals ',info
    WRITE(*,*)
    WRITE(*,"(3X,A)") 'EIGENVALUES:'
    WRITE(*,"(4F18.12)") eigval(:)
    WRITE(*,*)
    DEALLOCATE(eigval)
    DEALLOCATE(work)
    DEALLOCATE(matrix)
END PROGRAM eigenvalues