
MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, nwater, corr_time, nattype
  INTEGER                            :: nframes, nequil, nhist
  INTEGER                            :: nlayers, nskip, stride
  INTEGER                            :: maxneighbor, ind_rdf(2)
  INTEGER, ALLOCATABLE               :: moltype(:), atype(:), atype2(:)
  INTEGER, ALLOCATABLE               :: gofr(:) 
  INTEGER, ALLOCATABLE               :: coarse(:,:) 
  INTEGER, ALLOCATABLE               :: OH_index(:)
  INTEGER, ALLOCATABLE               :: ind_atom(:)
  INTEGER, ALLOCATABLE               :: neighborlist(:,:)
  REAL*8                             :: box(3), maxr, number_density
  REAL*8, ALLOCATABLE                :: pos(:,:), vel(:,:)
  REAL*8, ALLOCATABLE                :: layers(:) 
  REAL*8                             :: dt
  REAL*8, PARAMETER                  :: cutoff_OH = 1.3d0
  REAL*8, PARAMETER                  :: Pi=3.1415
END MODULE parameters

#ifdef __CUDA
MODULE cuda_parameters
  IMPLICIT NONE
  REAL*8, ALLOCATABLE, DEVICE        :: pos_d(:,:)
  REAL*8, DEVICE                     :: box_d(3)
  INTEGER, ALLOCATABLE, DEVICE       :: atype_d(:)
END MODULE cuda_parameters
#endif
  
SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE REMOVE_EQUIL_IPI
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+2)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL_IPI

SUBROUTINE REMOVE_EQUIL_VEL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+9)
    READ(2,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL_VEL

SUBROUTINE READ_ATOM 
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
  REAL*8                     :: box_tmp(2), box_min(3)
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
    box_min(iat)=box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, atype(ind), pos(1,ind), pos(2,ind), pos(3,ind)
  END DO

END SUBROUTINE READ_ATOM

SUBROUTINE READ_VELOCITIES 
  !Read LAMMPS custom velocities
  USE parameters, ONLY : vel, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
 
  DO iat=1,9
    READ(2,*)
  END DO

  DO iat = 1, natoms
    !READ(2,*) ind, atype(ind), vel(ind,1), vel(ind,2), vel(ind,3)
    READ(2,*) ind, atp, vel(1,ind), vel(2,ind), vel(3,ind)
  END DO

END SUBROUTINE READ_VELOCITIES

SUBROUTINE READ_ATOM_REDUCED 
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
  REAL*8                     :: box_tmp(2), box_min(3)
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
    box_min(iat)=box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, atype(ind), pos(1,ind), pos(2,ind), pos(3,ind)
    pos(:,ind) = pos(:,ind) * box + box_min 
  END DO

END SUBROUTINE READ_ATOM_REDUCED

SUBROUTINE READ_XYZ
  !Read xyz file. Line after natoms should contain box(1) box(2) box(3)
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  INTEGER                    :: get_ind_atom_from_xyz
  CHARACTER(5)               :: indatom
 
  READ(1,*)
  READ(1,*) box

  DO iat = 1, natoms
    READ(1,*) indatom, pos(1,iat), pos(2,iat), pos(3,iat)
    atype(iat)=get_ind_atom_from_xyz(TRIM(indatom))
  END DO

END SUBROUTINE READ_XYZ

SUBROUTINE READ_XYZ_IPI
  !Read xyz file from IPI
  USE parameters, ONLY : pos, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  CHARACTER(10)               :: trash
  DOUBLE PRECISION           :: tmp(3)
 
  READ(1,*)
  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) trash, pos(1,iat), pos(2,iat), pos(3,iat)
  END DO

END SUBROUTINE READ_XYZ_IPI

SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO
  USE parameters, ONLY : natoms, pos, box, atype
  IMPLICIT NONE
  INTEGER                              :: iat, ipol
  REAL*8                               :: mass
  REAL*8                               :: mtot, m, cm(3)
 
  cm=0.d0
  mtot=0.d0
 
  DO iat = 1,natoms

    m = mass(atype(iat))
    
    DO ipol = 1,3
      pos(ipol,iat) = pos(ipol,iat) & 
                    - nint(pos(ipol,iat)/box(ipol))*box(ipol)
    END DO  

    cm = cm + m*pos(:,iat)
    mtot = mtot + m

  END DO

  cm = cm / mtot
  DO iat = 1,natoms
    pos(:,iat) = pos(:,iat) - cm
  END DO

END SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO

SUBROUTINE COARSE_GRAIN_POS (frame)
  USE parameters, ONLY : natoms,pos, coarse, nlayers!, box, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  !DOUBLE PRECISION                     :: z
  INTEGER                              :: iat, ind_z, iwat
  INTEGER                              :: layer_index
  LOGICAL                              :: is_water_oxygen
  
  IF (nlayers < 2) RETURN

  coarse(frame,:) = 0
  iwat = 0

  DO iat = 1, natoms

    !Selecting only OW
    IF ( is_water_oxygen(iat) ) THEN

      iwat = iwat + 1
      ind_z = layer_index(pos(3,iat))
      coarse(frame,iwat) = ind_z

    END IF
    !IF (atype(iat)==1) then
    !  z = pos(3,iat) - box(3)/6.
    !  !z = z - nint(z/box(3))*box(3)
    !  !z = z - box(3)
    !  print *, z
    !END IF

  END DO

END SUBROUTINE COARSE_GRAIN_POS

SUBROUTINE READ_INDEX (index_file)
  USE parameters, ONLY : natoms, nframes, nequil, dt, &
                         ind_atom, nlayers, layers, nwater

  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file
  LOGICAL                    :: is_water_oxygen

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 3 lines of index file 
  READ(1, *) natoms, nframes, nequil, nlayers, dt

  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF

  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)") ind_atom(i)  
  END DO
 
  !Count number of water molecules and oxygen atoms
  nwater = 0
  DO i = 1, natoms
    IF (is_water_oxygen(i)) nwater = nwater+1
  END DO 

  CLOSE(1)

  nframes = nframes - nequil

END SUBROUTINE READ_INDEX

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jump between layers in the coarse grain function
  USE parameters, ONLY : coarse, nframes, nlayers
  IMPLICIT NONE
  INTEGER                    :: iw, frame, layer, tout, np
  INTEGER, PARAMETER         :: forgiv=20
  
  IF (nlayers < 2) RETURN

  np = size(coarse,2)

  DO iw = 1, np
  
    frame = 1
    tout  = 0
    layer = coarse(frame,iw)

    DO WHILE ( frame <= nframes ) 

      IF ( coarse(frame,iw) /= layer ) THEN
        tout = tout + 1
      ELSEIF ( tout > 0 ) THEN 
        coarse(frame-tout:frame,iw) = layer
        tout = 0
      END IF
  
      !If the particle really changed layer
      IF ( tout > forgiv ) THEN

        layer = coarse(frame,iw)
        tout = 0

      END IF

      frame = frame + 1

    END DO

  END DO

END SUBROUTINE SMOOTH_COARSE_GRAIN

SUBROUTINE APPLY_PBC(x_in,x_out)
  USE parameters, ONLY : box 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)        :: x_in(3)
  DOUBLE PRECISION                    :: x_out(3)
  INTEGER                             :: ipol

  DO ipol = 1,3

    x_out(ipol) = x_in(ipol) - nint(x_in(ipol)/box(ipol))*box(ipol) 
    x_out(ipol) = x_out(ipol) + box(ipol)/2.d0

  END DO

END SUBROUTINE APPLY_PBC 

SUBROUTINE DISTANCE_VECTOR( coord, wfc, d ) 
  USE parameters, ONLY : box
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: coord(3), wfc(3)
  DOUBLE PRECISION             :: d(4)
  INTEGER                      :: ipol
  
  DO ipol = 1,3
    d(ipol) = wfc(ipol) - coord(ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol))*box(ipol)
  END DO

  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR

SUBROUTINE IDENTIFY_OH_GROUPS 
  !For each H in the system, assign the index of the O atom it is bound to.
  !Return the indexes in the array OH_index
  USE parameters, ONLY : natoms, OH_index, box!, pos, atype
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2
  DOUBLE PRECISION           :: Dist, d, d_OH
  DOUBLE PRECISION           :: maxdist!, z
  LOGICAL                    :: is_hydrogen, is_oxygen
  INTEGER                    :: nH(natoms), ind_bkp(natoms), ind_out
  
  OH_index = 0
  nH=0

  DO iat1 = 1, natoms

    IF ( is_hydrogen(iat1) ) THEN

      d = 100.d0

      DO iat2 = 1, natoms

        IF ( is_oxygen(iat2) ) THEN

          d_OH = Dist(iat1,iat2)

          IF ( d_OH < d ) THEN

            ind_bkp(iat1) = OH_index(iat1)
            OH_index(iat1) = iat2 
            d = d_OH

          END IF

        END IF

      END DO

      nH(OH_index(iat1)) = nH(OH_index(iat1)) + 1

      IF ( d > 2 ) THEN
        print *, 'Hydrogen ', iat1, ' is not bound' 
      END IF

    END IF

  END DO

  !Check if a oxygen has more than 2 hydrogens
  DO iat1 = 1,natoms

    IF ( nH(iat1) > 2 ) THEN
  
      DO iat2 = 1, natoms

        maxdist=0
        IF ( OH_index(iat2) .eq. iat1 ) THEN
          IF ( Dist(iat1,iat2) > maxdist ) THEN
            ind_out = iat2
            maxdist = Dist(iat1,iat2)
          END IF
        END IF

      END DO

      !Assign 2nd nearest neighbor to the H farthest from oxygen iat1
      OH_index(ind_out) = ind_bkp(ind_out)
      
      IF ( nH(OH_index(ind_out)) .eq. 2 ) print *, "No way for index ", ind_out

    END IF

  END DO

  !DO iat1=1,natoms
  !  IF (atype(iat1)==1) then
  !    z = pos(3,iat1) - box(3)/6.
  !    z = z - nint(z/box(3))*box(3)
  !    print *, z
  !  END IF
  !END DO

END SUBROUTINE IDENTIFY_OH_GROUPS

SUBROUTINE DIST_UNIT_VECTOR(ind1,ind2,xyz,normvec)
  ! Unit vector from ind1 to ind2 
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  DOUBLE PRECISION                     :: xyz(3), normvec
  INTEGER                              :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(i,ind2) - pos(i,ind1)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  normvec=SQRT( SUM(xyz*xyz) )
  xyz = xyz / normvec 

END SUBROUTINE DIST_UNIT_VECTOR

SUBROUTINE MAKE_NEIGHBOR_LIST(atomtype1,atomtype2,rcut)
  !atomtypei=-1 means all atom types
  USE parameters, ONLY : natoms, atype, neighborlist, maxneighbor
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: atomtype1,atomtype2
  REAL*8, INTENT(IN)     :: rcut
  REAL*8                 :: Dist
  INTEGER                :: iat, jat, ineigh

  DO iat=1,natoms

    ineigh=0
    IF ((atype(iat)==atomtype1).or.(atomtype1.eq.-1)) THEN
      
      DO jat=1,natoms
   
        IF((atype(jat)==atomtype2).or.(atomtype2.eq.-1)) THEN

          IF (Dist(iat,jat)<rcut) THEN
            ineigh=ineigh+1
            neighborlist(iat,ineigh) = jat 
          END IF

        END IF
  
        IF (ineigh.eq.maxneighbor) EXIT

      END DO

    END IF

  END DO

END SUBROUTINE MAKE_NEIGHBOR_LIST

SUBROUTINE CROSSPROD( x, y, prod )
  !Cross product of 3D vectors x and y
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: prod
  REAL*8, INTENT(IN)         :: x(3),y(3)
  
  prod(1) = x(2)*y(3) - x(3)*y(2)
  prod(2) = x(3)*y(1) - x(1)*y(3)
  prod(3) = x(1)*y(2) - x(2)*y(1)

END SUBROUTINE CROSSPROD

SUBROUTINE Bubble_Sort(a,ind,n,ntot)
  !Sort array a in ascending order. Ind is an index array
  !a will contain the smallest n entries sorted from 1 to n
  INTEGER, INTENT(in) :: n, ntot
  REAL*8, DIMENSION(ntot) :: a
  INTEGER, DIMENSION(ntot) :: ind
  REAL*8 :: temp
  INTEGER :: i, j, itemp
  LOGICAL :: swapped
 
  !DO j = SIZE(a)-1, 1, -1
  DO j = 1, n 
    swapped = .FALSE.
    !DO i = 1, j
    DO i = ntot-1, j, -1
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        itemp = ind(i)
        ind(i) = ind(i+1)
        ind(i+1) = itemp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
RECURSIVE SUBROUTINE quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
END SUBROUTINE quicksort

REAL*8 FUNCTION mass(ind_atom)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                   :: ind_atom
  
  SELECT CASE (ind_atom)
    CASE (1)
      mass = 47.867
    CASE (2)
      mass = 2.018
    CASE (3)
      mass = 15.999
    CASE (4) 
      mass = 15.999
  END SELECT

END FUNCTION mass

LOGICAL FUNCTION is_water_oxygen(iat)
  USE parameters, only: atype!,ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: iat

  is_water_oxygen=.false.
  !IF (ind_atom(iat).eq.4) is_water_oxygen=.true.
  IF (atype(iat).eq.1) is_water_oxygen=.true.

END FUNCTION is_water_oxygen

INTEGER FUNCTION layer_index(pos_z)
  !For layer resolved analysis - Get the index of Ow with z coordinate
  USE parameters, ONLY : box, nlayers, layers
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: pos_z
  DOUBLE PRECISION                     :: ind, z
  INTEGER                              :: i
 
  IF (nlayers < 2) THEN
    layer_index = 1
    RETURN
  END IF

  z = pos_z - box(3)/6.
  z = z - nint( z / box(3) ) * box(3)
  !print *,z
  !z = z + box(3)/2.d0
  i = 0 

  DO WHILE ( i < nlayers ) 
 
    i = i + 1
    IF ( ( z > layers(i) ) .and. ( z <= layers(i+1) ) ) THEN 
      layer_index = i
      RETURN
    END IF

  END DO

  layer_index = 0 
  PRINT *, "Could not find layer for particle with z ", z
  STOP

END FUNCTION layer_index

LOGICAL FUNCTION in_layer( iat, ti, tf, sep )
  !Test if atom is inside a layer from time ti to tf
  USE parameters, ONLY : coarse, nlayers
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: iat, ti, tf
  INTEGER                    :: layer, ts
  INTEGER, OPTIONAL          :: sep

  ts = 1
  IF (present(sep)) ts=sep

  IF ( nlayers > 1 ) THEN

    layer = coarse(ti,iat)
    in_layer = ALL( coarse(ti:tf:ts,iat) == layer )

  ELSE 

    in_layer = .true.

  END IF

END FUNCTION in_layer

DOUBLE PRECISION FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  DOUBLE PRECISION                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(i,ind1) - pos(i,ind2)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

DOUBLE PRECISION FUNCTION Dist_Center(center,coord)
  ! Distance from a point coord to the center
  USE parameters, ONLY :  box
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: center(3),coord(3)
  DOUBLE PRECISION                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  xyz=0.0
  DO i = 1,2

    xyz(i) = coord(i) - center(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist_Center = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist_Center

SUBROUTINE Vector_Distance(a,b,d_out)
  !Vector distance between vectors a and b. 
  USE parameters, only : box
  IMPLICIT NONE
  REAL*8, INTENT(in) :: a(3), b(3)
  REAL*8 :: d_out(4)
  INTEGER :: pol

  DO pol=1,3
    d_out(pol) = a(pol)-b(pol)
    d_out(pol) = d_out(pol) - nint(d_out(pol)/box(pol))*box(pol)
  END DO

  d_out(4) = sqrt(sum(d_out(1:3)**2))

END SUBROUTINE Vector_Distance

SUBROUTINE Vector_Distancei(a,b,d_out)
  !Vector distance between pos vectors with indexes a and b. 
  USE parameters, only : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(in) :: a, b
  REAL*8 :: d_out(4)
  INTEGER :: pol

  CALL Vector_Distance(pos(:,a),pos(:,b),d_out)

END SUBROUTINE Vector_Distancei

LOGICAL FUNCTION is_oxygen(ind)
  USE parameters, ONLY : atype !ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  !IF ( ( ind_atom(ind) == 3 ) .or. ( ind_atom(ind) == 4 ) ) THEN
  IF  ( atype(ind) == 1 ) THEN
    is_oxygen = .true.
  ELSE
    is_oxygen = .false.
  END IF

END FUNCTION is_oxygen

LOGICAL FUNCTION is_hydrogen(ind)
  USE parameters, ONLY : atype !ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  !IF ( ind_atom(ind) == 2 ) THEN
  IF ( atype(ind) == 2 ) THEN
    is_hydrogen = .true.
  ELSE
    is_hydrogen = .false.
  END IF

END FUNCTION is_hydrogen

REAL*8 FUNCTION Angle(C, a1, a2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  IMPLICIT NONE
  REAL*8, INTENT(IN), DIMENSION(3)       :: C, a1, a2
  REAL*8 :: d1(4), d2(4)

  CALL Vector_Distance(C,a1,d1)
  CALL Vector_Distance(C,a2,d2)
  Angle = dot_product(d1(1:3),d2(1:3)) / ( d1(4) * d2(4) )

END FUNCTION Angle

REAL*8 FUNCTION Angle_reference(C, a1, reference)
  !Angle between C-a1 and reference
  !reference should be a unit vector
  !Return value in cos(Angle)
  IMPLICIT NONE
  REAL*8, INTENT(IN), DIMENSION(3)       :: C, a1, reference
  REAL*8 :: d1(4)

  CALL Vector_Distance(C,a1,d1)
  Angle_reference = dot_product(d1(1:3),reference) / d1(4)

END FUNCTION Angle_reference

REAL*8 FUNCTION Anglei(iC, ia1, ia2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters, only : pos
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iC, ia1, ia2
  REAL*8 :: Angle

  Anglei = Angle(pos(:,iC),pos(:,ia1),pos(:,ia2)) 

END FUNCTION Anglei

REAL*8 FUNCTION Angle_Bissector(C, a1, reference)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters, only : pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)       :: C, a1(2)
  REAL*8 :: reference(3)
  REAL*8 :: d1(4), d2(4)

  CALL Vector_Distancei(C,a1(1),d1)
  CALL Vector_Distancei(C,a1(2),d2)
  d1 = (d1+d2)/2.
  d1(4) = sqrt(sum(d1(1:3)*d1(1:3)))
  !CALL Vector_Distance(pos(:,C),reference,d2)
  !Angle_Bissector = dot_product(d1(1:3),d2(1:3)) / ( d1(4) * d2(4) )
  Angle_Bissector = dot_product(d1(1:3),reference) / d1(4) 

END FUNCTION Angle_Bissector

INTEGER FUNCTION Coordination_Number(iat,neighbor_type,cutoff)
  !This function assumes that a neighborlist was created
  USE parameters, ONLY: natoms, atype, neighborlist, maxneighbor
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, neighbor_type 
  REAL*8, INTENT(IN)         :: cutoff
  INTEGER                    :: i, cn
  REAL*8                     :: Dist

  cn=0
  DO i=1,maxneighbor
    if (neighborlist(iat,i).eq.0) EXIT
    if (atype(neighborlist(iat,i)).eq.neighbor_type) then
      if (Dist(iat,neighborlist(iat,i))<cutoff) cn=cn+1
    end if
  END DO

  Coordination_Number=cn

END FUNCTION Coordination_Number

INTEGER FUNCTION get_ind_atom_from_xyz(indatom)
  IMPLICIT NONE
  CHARACTER(5), INTENT(IN) :: indatom
  INTEGER                  :: ind

  SELECT CASE (indatom)
    CASE("Ti") 
      ind=1
    CASE("H")  
      ind=2
    CASE("O")  
      ind=3
    CASE("Na") 
      ind=4
    CASE("Cl") 
      ind=5
    CASE DEFAULT 
      ind=0
  END SELECT

  get_ind_atom_from_xyz=ind

END FUNCTION get_ind_atom_from_xyz
