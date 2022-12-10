PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  CALL DEFINE_ADDITIONAL_ATOM_TYPES
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) CALL COMPUTE_RDF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, number_density, &
                         nhist, gofr, nframes, nequil, &
                         box, ind_rdf, maxr, atype, stride
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist
  READ(1,*) box !3 dimensions
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species

  maxr=minval(box)/2. !Max distance considered in g(r)
  number_density = 0.d0

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(gofr(nhist)); gofr=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE DEFINE_ADDITIONAL_ATOM_TYPES
  USE parameters, ONLY : natoms, atype, atype2
  IMPLICIT NONE
  INTEGER :: iat, coordination_number_, CN
  REAL*8, PARAMETER :: rcut=2.2 !Angstroms

  CALL READ_ATOM_REDUCED
  ALLOCATE(atype2(natoms))
  atype2=atype

  DO iat=1,natoms
    IF (atype(iat).eq.3) THEN !Define O atoms first
      CN = coordination_number_(iat,1,rcut)
      IF (CN.eq.2) THEN
        atype2(iat) = 32
        atype(iat)=30
      ELSE IF (CN.eq.3) THEN
        atype2(iat) = 33
        atype(iat)=30
      ELSE IF (CN.lt.2) THEN
        atype2(iat) = 31
      END IF
    END IF
  END DO

  DO iat=1,natoms
    IF (atype(iat).eq.1) THEN !Define Ti atoms
      CN = coordination_number_(iat,30,rcut)
      IF (CN.eq.5) atype2(iat)=15
      IF (CN.eq.6) atype2(iat)=16
    END IF
  END DO

  REWIND(1)

END SUBROUTINE DEFINE_ADDITIONAL_ATOM_TYPES

SUBROUTINE COMPUTE_RDF
  USE parameters, ONLY : natoms, pos, ind_rdf, box, &
                         nhist, gofr, atype2, maxr, number_density
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: n_coord 
  INTEGER                     :: ind, gofr_tmp(nhist), n_tmp
  INTEGER                     :: n
  REAL*8                      :: Dist, d

  n = 0 !total number of pairs
  DO iat = 1, natoms

    IF ( ( atype2(iat) == ind_rdf(1) ) .or. ( ind_rdf(1) == 0 ) ) THEN

      gofr_tmp=0
      n_tmp=0
!$omp parallel do private(d,ind) reduction(+:gofr_tmp,n_tmp)
      DO iat2 = 1, natoms

        IF ( ( ( atype2(iat2) == ind_rdf(2) ) .or. ( ind_rdf(2) == 0 ) ) &
             .and. ( iat.ne.iat2 ) ) THEN

          d = Dist(iat,iat2)
          ind = int(nhist*d/maxr) + 1
          IF (ind<=nhist) gofr_tmp(ind) = gofr_tmp(ind) + 1
          n_tmp = n_tmp + 1

        END IF

      END DO
!$omp end parallel do
      gofr = gofr + gofr_tmp
      n = n + n_tmp 

    END IF

  END DO

  number_density = number_density + float(n)/product(box)

END SUBROUTINE COMPUTE_RDF

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, gofr, maxr, number_density, &
                         atype2, natoms, ind_rdf, nframes
  IMPLICIT NONE
  REAL*8                     :: bin, aver, vol
  REAL*8, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist, N 
  INTEGER                    :: get_natoms_type
  
  OPEN(unit = 2,file = "rdf.dat")

  write(2,*) '#', 'Z', 'RDF'

  N = get_natoms_type(atype2,natoms,ind_rdf(1)) 

  bin =  maxr/float(nhist)
  !number_density = number_density !/ float(nframes)
  number_density = float(nframes) * float(N)


  DO ihist = 1, nhist

    vol = FourPioverThree * ( (float(ihist)*bin)**3 - ((float(ihist)-1)*bin)**3 )
    aver = float(gofr(ihist)) / number_density / vol
    WRITE(2,*) bin*(float(ihist)-0.5), aver 

  END DO

  CLOSE(2)

END SUBROUTINE

INTEGER FUNCTION get_natoms_type(atype,natoms,ind_rdf)
  !Return natoms with type ind_rdf
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: natoms, atype(natoms), ind_rdf
  INTEGER                    :: iat, irdf, N

  N = 0

  DO iat = 1, natoms

    IF ( ( atype(iat).eq.ind_rdf ) .or. &
       ( ind_rdf.eq.0 ) ) N = N + 1 

  END DO

  get_natoms_type = N

END FUNCTION get_natoms_type

INTEGER FUNCTION coordination_number_(iat,atyp,rcut)
  USE parameters, only : natoms,atype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iat,atyp
  REAL*8, INTENT(IN) :: rcut
  INTEGER :: i, cn
  REAL*8 :: d, Dist

  cn = 0 
  DO i=1,natoms

    IF ((atype(i).eq.atyp).and.(i.ne.iat)) THEN
      d = Dist(i,iat)
      IF (d<rcut) cn = cn + 1
    END IF

  END DO

  coordination_number_=cn

END FUNCTION coordination_number_
