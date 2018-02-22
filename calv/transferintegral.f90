
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

  integer ii,jj,nb,nb1,nb2,imol1,imol2,info,norbt(2),i1,i2
  integer,allocatable :: ipiv(:)
  character*20 flog1,flog2,dimerlog,mol1,mol2
  double precision c1fc2,DDOT,el1,el2,eh1,eh2
  double precision Jll, Jlh,Jhl,Jhh,sl1,sl2,sh1,sh2,sll,slh,shl,shh
  double precision  tll,thh,tlh,thl
  double precision,allocatable :: f(:,:),S(:,:),e(:),c(:,:),cl1(:),cl2(:),ctmp(:),work(:),SS(:,:),ch1(:),ch2(:)
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
  allocate(ch1(nb),cl1(nb),ch2(nb),cl2(nb))  ! add by Yecheng Zhou

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
  allocate(c(nb1,nb1))
  allocate(e(nb1))

  do while(.True.)
    read(99,'(a11)') flagc
    if(flagc.eq.' Alpha MOs:') exit
  enddo
  call read_MO1(e,c,nb1,99)
  close(99)
!write(*,*) "read mol 2  finished"


  !%%%%%%%%%%!
  !% cal c1 %!
  !%%%%%%%%%%!

  cl1=0.0d0; 
  ch1=0.0d0;  !! add by Yecheng Zhou
  do ii=1,nb1
    ch1(ii)=c(ii,imol1)
    cl1(ii)=c(ii,imol1+1) !! add by Yecheng Zhou
  enddo
  deallocate(e,c)



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
  allocate(c(nb2,nb2))
  allocate(e(nb2))

  do while(.True.)
    read(99,'(a11)') flagc
    if(flagc.eq.' Alpha MOs:') exit
  enddo
  call read_MO1(e,c,nb2,99)
  close(99)

!write(*,*) "read mol 1-2  finished"

  !%%%%%%%%%%!
  !% cal c2 %!
  !%%%%%%%%%%!

  cl2=0.0d0
   ch2=0.0d0
  do ii=1,nb2
    ch2(ii+nb1)=c(ii,imol2)
    cl2(ii+nb1)=c(ii,imol2+1)
  enddo
  deallocate(e,c)



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% calculate transfer integral %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  el1=0.0d0
  eh1=0.0d0
  do ii=1,nb
    do jj=1,nb
      el1=el1+cl1(ii)*F(ii,jj)*cl1(jj)
      eh1=eh1+ch1(ii)*F(ii,jj)*ch1(jj)
    enddo
  enddo
  
  el2=0.0d0
  eh2=0.0d0
  do ii=1,nb
    do jj=1,nb
      el2=el2+cl2(ii)*F(ii,jj)*cl2(jj)
      eh2=eh2+ch2(ii)*F(ii,jj)*ch2(jj)
    enddo
  enddo

  Jll=0.0d0
  Jlh=0.0d0
  Jhl=0.0d0
  Jhh=0.0d0
  do ii=1,nb
    do jj=1,nb
      Jll=Jll+cl1(ii)*F(ii,jj)*cl2(jj)
      Jlh=Jlh+cl1(ii)*F(ii,jj)*ch2(jj)
      Jhl=Jhl+ch1(ii)*F(ii,jj)*cl2(jj)
      Jhh=Jhh+ch1(ii)*F(ii,jj)*ch2(jj)
    enddo
  enddo

  sl1=0.0d0
  sh1=0.0d0
  do ii=1,nb
    do jj=1,nb
      sl1=sl1+cl1(ii)*S(ii,jj)*cl1(jj)
      sh1=sh1+ch1(ii)*S(ii,jj)*ch1(jj)
    enddo
  enddo

  sl2=0.0d0
  sh2=0.0d0
  do ii=1,nb
    do jj=1,nb
      sl2=sl2+cl2(ii)*S(ii,jj)*cl2(jj)
      sh2=sh2+ch2(ii)*S(ii,jj)*ch2(jj)
    enddo
  enddo

  sll=0.0d0
  slh=0.0d0
  shl=0.0d0
  shh=0.0d0
  do ii=1,nb
    do jj=1,nb
      Sll=Sll+cl1(ii)*S(ii,jj)*cl2(jj)
      Slh=Slh+cl1(ii)*S(ii,jj)*ch2(jj)
      Shl=Shl+ch1(ii)*S(ii,jj)*cl2(jj)
      Shh=Shh+ch1(ii)*S(ii,jj)*ch2(jj)
    enddo
  enddo

  tll=(Jll-0.5d0*(el1+el2)*Sll)/(1-Sll*Sll)*27.2116d0*1000d0
  thh=(Jhh-0.5d0*(eh1+eh2)*Shh)/(1-Shh*Shh)*27.2116d0*1000d0
  thl=(Jhl-0.5d0*(eh1+el2)*Shl)/(1-Shl*Shl)*27.2116d0*1000d0
  tlh=(Jlh-0.5d0*(el1+eh2)*Slh)/(1-Slh*Slh)*27.2116d0*1000d0 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% save transferintegral to OUTPUT %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!  write(6,'(2(xi0),f12.5)') imol1,imol2, t 
  write(6,'(2(xa),4(xf12.5))') "Tll Thh Thl Tlh for ", trim(dimerlog),Tll, Thh, Thl, Tlh 
!#  close(10)
!end do
!  close(10)
!  9994 write(6,'(a)') 'All of the integral calculation finished'

  !%%%%%%%%%%%%%%!
  !% normal end %!
  !%%%%%%%%%%%%%%!

  deallocate(cl1,cl2,ctmp,F,S,ch1,ch2)
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

