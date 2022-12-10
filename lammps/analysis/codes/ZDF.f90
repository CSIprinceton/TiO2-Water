!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: nhist
  INTEGER, ALLOCATABLE                     :: hist(:,:)
  INTEGER                                  :: ntype, zdir
  LOGICAL                                  :: input_is_xyz=.false.
  REAL*8                                   :: zoffset
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  USE histogram, only: input_is_xyz
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL(input_is_xyz)
  
  DO frame = 1,nframes
    if (input_is_xyz) then
      CALL READ_XYZ
    else
      CALL READ_ATOM_REDUCED
    end if
    IF (MOD(frame,nskip)==0) THEN
      !CALL SET_CENTER_OF_MASS_TO_ZERO
      CALL MAKE_HISTOGRAM
    END IF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, ntype, nhist, input_is_xyz, zdir, zoffset
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, nskip, nhist, zdir
  READ(1,*) zoffset


  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist(ntype,nhist)); hist=0

  !Read the index file
  IF (input_is_xyz) THEN
    READ(1,*) box 
    DO i = 1,natoms
      READ(1,fmt = "(I2)") atype(i)  
    END DO
  END IF

  CLOSE(1)
  OPEN(unit  =  1,file  =  pos_file)
 
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram, ONLY  : hist, nhist, zdir, zoffset
  USE parameters, ONLY : pos, box, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: i, ind
  REAL*8           :: z

  DO i = 1,natoms

      !Applying PBC
      z = pos(zdir,i) + zoffset
      z = z  - nint(z/box(zdir))*box(zdir) 
      z = z + box(zdir)/2. ! for data analysis only

      !Assigning the index to the histogram
      ind = int( z * float(nhist) / box(zdir) ) + 1
      hist(atype(i),ind) = hist(atype(i),ind) + 1

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nhist, hist, zdir
  USE parameters, ONLY : box, nframes, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: bin,dV

  OPEN(unit = 2,file = "ZDF.dat")

  bin = box(zdir)/float(nhist)
  dV = product(box)/box(zdir)*bin

    DO i = 1,nhist
      WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + float(i-1)*bin - box(zdir)/2.d0, &
                           & float(hist(:,i))/(0.03316d0*dV*float(nframes)/float(nskip))
      !WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + float(i-1)*bin, &
      !                     & float(hist(:,i))/(0.03316d0*dV*float(nframes/nskip))
    END DO
    !DO i = 1,nhist
    !  WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + float(nhist+i-1)*bin - box(3)/2.d0, &
    !                       & float(hist(:,i))/(0.03316d0*dV*float(nframes))
    !END DO

  CLOSE(2)

END SUBROUTINE
