! ifort overlap.f90  -mkl
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% calculate transfer integral with the onsite energy correction method %!
!% according to Bredas,JACS,2006-128-9882                               %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%          J12-0.5*(e1+e2)*S12                                         %!
!% J12_eff=---------------------                                        %!
!%                 1-S12^2                                              %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% ljwang,Feb.16,2009     modified and expanded by Yecheng Zhou         %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

program pr_transferintegral
  implicit none

  !%%%%%%%%%%%%%%!
  !% parameters %!
  !%%%%%%%%%%%%%%!

  integer ii,jj,nb,nb1,nb2,imol1,imol2,info,norbt(2),i1,i2,orbt1,orbt2,nmax
  integer,allocatable :: ipiv(:)
  double precision c1fc2,DDOT,smax,overlap
  double precision,allocatable :: f(:,:),S(:,:),e(:),c(:,:),ctmp(:),work(:),SS(:,:)
  double precision,allocatable :: c1all(:,:),c2all(:,:),c1(:),c2(:) ! for HOMO-10 to HOMO and LUMO to lumo+10,
  character*20 flog1,flog2,dimerlog,mol1,mol2  ! HOMO of these molecules
  character*100 flagc 
  external DGEMM,DGETRF,DGETRI,DGEMV,DSCAL,DDOT


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% get the name of dimerlog,flog1,flog2 %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
  if (iargc()/=5) then 
     write(*,*) "input parameter missing,stop"!    ch1(ii)=c(ii,imol1)
     write(*,*) "first-log, second-log, dimer-log, first-HOMO, second-HOMO" ! cl1(ii)=c(ii,imol1+1) !! add by Yecheng Zhou
     stop
  end if 
  
  call getarg(1,flog1)
  call getarg(2,flog2) 
  call getarg(3,dimerlog)  
  call getarg(4,mol1) 
  call getarg(5,mol2)
  read(mol1,"(i20)") i1 
  read(mol2,"(i20)") i2 
  imol1=i1
  imol2=i2 
 
  if (i1<10 .and. i2 <10) then  !i1 and i2 are the molecule type rather than orbital numebr 
      open(100,file="in",status='old') 
      read(100,*) !flname 
      read(100,*) !ncpu 
      read(100,*) !cutoff 
      read(100,*) (norbt(ii),ii=1,2)
      close(100)
      imol1= norbt(i1)
      imol2= norbt(i2)
   end if 
!  write(*,*) trim(flog1),"x",trim(flog2),"x",trim(dimerlog),trim(mol1),trim(mol2)," x ",imol1,imol2
!  write(*,*) flog1, flog2, dimerlog 

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% read information from INPUT %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!#  read(5,*) imol1
!#  read(5,*) imol2 
     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% read information from dimerlog %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  open(unit=99,file=dimerlog,status='old',err=9993)

  do while(.True.)
    read(99,'(a12)') flagc
    if(flagc.eq.'    NBasis =') exit
  enddo
  backspace(99)
  read(99,'(A12,I4)') flagc,nb
!   write(*,*) "read Nbasis ok"
   
  allocate(f(nb,nb),s(nb,nb),c(nb,nb),e(nb),ipiv(nb),work(nb),ctmp(nb),SS(nb,nb)) ! deleted c1 and c2 replace by following 
  allocate(c1all(nb,20),c2all(nb,20))  ! add by Yecheng Zhou, 10 HOMO and LUMO

  do while(.True.)
    read(99,'(a16)') flagc
    if(flagc.eq.' *** Overlap ***') exit
  enddo
  call read_lowtri(s,nb,99)
!  write(*,*) "read overlap ok"

  do while(.True.)
    read(99,'(a11)') flagc
    if(flagc.eq.' Alpha MOs:') exit
  enddo
  call read_MO1(e,c,nb,99)
!  write(*,*) "read Alpha MOs: ok"
  close(99)

!write(*,*) "read mol 1  finished"

  !%%%%%%%%%%%%%%%%%%%!
  !% cal Fock matrix %!
  !%%%%%%%%%%%%%%%%%%%!

  ss=s
  call DGEMM('N','N',nb,nb,nb,1.0d0,ss,nb,c,nb,0.0d0,f,nb)
  do ii=1,nb
    call DSCAL(nb,e(ii),f(1,ii),1)
  end do

  call DGETRF(nb,nb,c,nb,ipiv,info)
  call DGETRI(nb,c,nb,ipiv,work,nb,info)
  ss=f
  call DGEMM('N','N',nb,nb,nb,1.0d0,ss,nb,c,nb,0.0d0,f,nb)
  deallocate(e,c,ipiv,work,ss)


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% read information from flog1 %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  open(unit=99,file=flog1,status='old',err=9991)
  
  do while(.True.)
    read(99,'(a12)') flagc
    if(flagc.eq.'    NBasis =') exit
  enddo
  backspace(99)
  read(99,"(A12,I4)") flagc,nb1
  allocate(c(nb1,nb1),c1(nb1))
  allocate(e(nb1))

  do while(.True.)
    read(99,'(a11)') flagc
    if(flagc.eq.' Alpha MOs:') exit
  enddo
  call read_MO1(e,c,nb1,99)
  close(99)
  c1all(1:nb1,1:2) = c(1:nb1,imol1-9:imol1+10)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% read information from flog2 %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  open(unit=99,file=flog2,status='old',err=9992)
  do while(.True.)
    read(99,'(a12)') flagc
    if(flagc.eq.'    NBasis =') exit
  enddo
  backspace(99)
  read(99,"(A12,I4)") flagc,nb2
  allocate(c(nb2,nb2),c1(nb1))
  allocate(e(nb2))

  do while(.True.)
    read(99,'(a11)') flagc
    if(flagc.eq.' Alpha MOs:') exit
  enddo
  call read_MO1(e,c,nb2,99)
  close(99)
  c2all(1:nb2,1:2) = c(1:nb2,imol1-9:imol1+10)

!####################
! calculate overlap for these orbitals
!###################
  do orbt1=1,20
    c1=c1all(:,orbt1)
    smax=0.0;nmax=orbt1
    do orbt2=1,20
      overlap=0.0
      c2=c1all(:,orbt2)
      do ii=1,nb
        do jj=1,nb
          overlap=overlap+c1(ii)*S(ii,jj)*c2(jj)
        enddo
      enddo
      if (overlap >=smax) then
        smax=overlap
        nmax=imol2-10+orbt2  ! 1 for imol2-9
      end if 
    end do 
    write(6,'(a,xi5,xi5,xf12.5)') "orital1 max overlap with orbital2 and its overlap", imol1-10+orbt2,nmax,smax
  end do

  deallocate(c1all,c2all,c1,c2,c,ctmp,F,S)
  stop

  !%%%%%%%%%%%%%!
  !% error end %!
  !%%%%%%%%%%%%%!

  9990 write(6,'(a)') '***** FILE input NOT FOUND! *****'
  stop
  9991 write(6,'(a)') '***** FILE'//flog1(1:len_trim(flog1))//' NOT FOUND! *****'
  stop
  9992 write(6,'(a)') '***** FILE'//flog2(1:len_trim(flog2))//' NOT FOUND! *****'
  stop
  9993 write(6,'(a)') '***** FILE'//dimerlog(1:len_trim(dimerlog))//' NOT FOUND! *****'
  stop

end program


!!###########################################################################################################
!%%%%%%%%%%%%%%%%!
!% read low tri %!
!%%%%%%%%%%%%%%%%!

subroutine read_lowtri(a,n,iunt)
  implicit none

  integer n,iunt,ii,jj,kk,nblk,npblk,jb,je,kb,ke,idum
  double precision a(n,n)

  nblk=5
  npblk=(n-1)/nblk+1
  do ii=1,npblk
    read(iunt,*)
    jb=(ii-1)*nblk+1
    je=n
    kb=jb
    do jj=jb,je
      ke=jb+nblk-1
      if(jj<ke) ke=jj
      read(iunt,*) idum,(a(jj,kk),kk=kb,ke)
      do kk=kb,ke
        a(kk,jj)=a(jj,kk)
      enddo
    enddo
  enddo

  return
end subroutine


!%%%%%%%%%%%!
!% read MO %!
!%%%%%%%%%%%!

subroutine read_MO(e,c,n,iunt)
  implicit none

  integer n,iunt,ii,jj,kk,jb,je,nblk,npblk
  double precision e(n),c(n,n)
  character*80 cha

  nblk=5
  npblk=(n-1)/nblk+1
  do ii=1,npblk
    read(iunt,*) cha
    read(iunt,*) cha

    jb=(ii-1)*nblk+1
    je=ii*nblk
    if(je>n) je=n
    read(iunt,'(A21,5f10.5)') cha,(e(jj),jj=jb,je)
    do kk=1,n
      read(iunt,'(A21,5f10.5)') cha,(c(kk,jj),jj=jb,je)
    enddo
  enddo

  return
end subroutine


!%%%%%%%%%%%%!
!% read MO1 %!
!%%%%%%%%%%%%!

subroutine read_MO1(e,c,n,iunt)
  implicit none

  integer n,iunt,ii,jj,kk,jb,je,nblk,npblk,nline
  double precision e(n),c(n,n)
  character*80 cha

  nblk=5
  npblk=(n-1)/nblk+1
  nline=0 
!  write(*,*) npblk,floor(real((n-1)/nblk))+1
  do ii=1,npblk
    read(iunt,*) cha
    nline=nline+1
    jb=(ii-1)*nblk+1
    je=ii*nblk
    if(je>n) je=n
    read(iunt,'(A21,5f10.5)',err=200) cha,(e(jj),jj=jb,je)
    nline=nline+1
    do kk=1,n
      read(iunt,'(A21,5f10.5)',err=210) cha,(c(kk,jj),jj=jb,je)
      nline=nline+1
    enddo
  enddo
   goto 220 
200  write(*,*) "read e error", jj 
210    backspace(iunt)
     write(*,*)  "read AMOs error", jj, kk,nline+1
     read(iunt,*) cha
     write(*,*) cha
     read(iunt,*) cha
     write(*,*) cha
220  return
end subroutine

