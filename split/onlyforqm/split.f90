program main 
!ccc molecule-atom-distance-basis-core cccccccccccccccccc
!use lib 
implicit none 
integer::nam(10000),nm,i,j,k,l,m,n,count(10000),ii,rt,ri,tm,iin,iim  !number of atoms in molecule, and numbers of molecules
real::c(3,10000,10000),dc(3,10000),cen(3,10000),dis(3),cutoff,dtc(3,10000) !c(im,ia,3)
real::tdis,cell(3,3),kct,tdis1,tdis2,vtmp(3),vtmp0(3),dis1(3),laxs(3),tv1(3),tdis3
character(2)::elemt(10000,10000),delemt(20000) ! ,ALLOCATABLE
character(2),allocatable::tpielemt(:,:,:),pielemt(:,:,:)
character(20)::flname,mnum,mnum2,flnm,itpfile
integer::km(10000),ncpu,tmp,tmp2,piindx(4,100),jj,mindx(4),piindxnew(4,100) !,pair(1000,1000,2)
character(20)::sets 
character(6)::label(1000)
character(200)::cmd
integer::nhomo,mhomo,rc(3),frc(3),nchain,chain(1000,10),ichain(10),tnv
logical::levelm(10000),stt(10000,10000)
real::tmp_pi(3,1000),sc
real,allocatable::picrd(:,:,:,:),tpicrd(:,:,:,:),cip(:,:),tpichg(:,:,:)
real::outparam(10),chg(1000),outparam1(10),outparam2(10),tchg(1000)
integer::anglhist(100),rhist(100),rchist(100)
integer::kkm,kkn,nnc,nnca(100)
character(15)::method 
logical::debug(5),oldcon

debug=.False.

open(100,file="in",status='old') 
read(100,*) flname
read(100,*) ncpu 
read(100,*) cutoff 
cutoff=cutoff*cutoff
read(100,*) 
read(100,"(a20)") sets
close(100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read molecules coordinate           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dc=0
call readgro(flname,cell,c,elemt,nm,nam,km) ! coordinate in nm
write(*,"(a,9(xf12.3))") "cell", cell
tv1(1)=cell(1,1)
tv1(2)=cell(2,2)
tv1(3)=cell(3,3)

do m=1,nm 
   do n=1,nam(m) 
     if (n/=1) then 
           rc=0
           do i=1,3
            if (c(i,n,m)-c(i,n-1,m)> 0.5*tv1(i)) rc(i)=-1
            if (c(i,n,m)-c(i,n-1,m)<-0.5*tv1(i)) rc(i)= 1
           end do !i
           do i=1,3
               c(1:3,n,m)=c(1:3,n,m)+rc(i)*cell(1:3,i)
           end do !i
     end if 
  end do ! n
end do ! m 

!##################################################
      !     calculate center coord of molecule    !
!##################################################
kct=0.d0;        cen=0.d0;
 write(*,*) "number of molecules: ", nm 
do m=1,nm 
   do n=1,nam(m) 
       cen(:,m)=cen(:,m)+c(:,n,m)
   end do 
   cen(:,m)=cen(:,m)/nam(m)
!    Write(*,*) m , "central calculation finished" 
end do 
!Write(*,*) "central calculation finished" 

if (debug(1)) stop 
!##################################################
!    select  dimer 
!##################################################
write(*,*) "select dimers "
inquire(file="connection.dat", exist=oldcon)

if (oldcon .eqv. .false.) then 
dis=10000 
open(131, file = "connection.dat", status = 'replace') 
!      write(131,*) nm ! number of molecules 
do m=1,nm 
    count(m)=0
    do n=m+1,nm
        tdis1=100
        do iim=1,nam(m)
        do iin=1,nam(n)
            vtmp(1:3)=c(1:3,iin,n)
            vtmp0(1:3)=c(1:3,iim,m)
            rc=0
            do i=1,3
             if (vtmp(i)-vtmp0(i) >  0.5*tv1(i)) rc(i)=-1
             if (vtmp(i)-vtmp0(i) < -0.5*tv1(i)) rc(i)= 1
            end do !i
            do i=1,3
                vtmp(1:3)=vtmp(1:3)+rc(i)*cell(1:3,i)
            end do !i
            vtmp(1:3)=vtmp(1:3)-vtmp0(1:3)
            tdis=vtmp(1)*vtmp(1)+vtmp(2)*vtmp(2)+vtmp(3)*vtmp(3)
            if (tdis<tdis1) then 
              tdis1=tdis  ! find minimum to tdis1 
              frc=rc
            end if 
        enddo !iin
        end do !iim
!##################################################
!    selection based on nearest atom distance.only calculate dimers with distance less than 5A
!##################################################
        if (tdis1 < cutoff) then 
            count(m)=count(m)+1
            vtmp(1:3)=cen(1:3,n)
            do i=1,3
              vtmp(1:3)=vtmp(1:3)+frc(i)*cell(1:3,i)
            end do !i
            vtmp(1:3)=vtmp(1:3)-cen(1:3,m)  
           write(131,"(i5,3i3,i5,4(xf10.4))") m, frc(1:3),n,(vtmp(i),i=1,3),sqrt(tdis1) 
        end if 
    end do !n
end do !m
close(131)

else !!!!! connection.dat exist
write(*,*) "connection.dat exist, no com and run file will be writen"
debug(2)=.True. 
count=0
end if

!##################################################
!    Generate rw input file
!##################################################
write(*,*) "reading connection.dat "
open(131, file = "connection.dat", status = 'old') 
open(191, file ='grp', status = 'replace') !gather the results from integral calculation 
write(191,"(a)") "#!/bin/bash "

open(211, file ='tcon', status = 'replace') !  information to run random walk simulation
write(211,'(i0,2xa,2xe15.5)') nm,'1000000 ',kct ! information to run random walk simulation 

rt=10000 ! initial time > 32, to create a new run scripts
ri=1
levelm=.true.
stt=.false.

write(*,*) "begin output g16 and descriptors"
do m=1,nm 
!##################################################
!    Generate list of job need to run 
!##################################################
  if (rt>32) then   ! used to split run file, count number of calculations
      write(flname,"(a,i0)") "run",ri 
      if (debug(2) .eqv. .false.) open(181, file =flname, status = 'replace') !command to run g16
      if (debug(2) .eqv. .false.) write(181,"(a)") "#!/bin/bash "
      rt=0; 
  end if 
  
  write(211,'(i0,xi0,xi0,3(xe13.5))') m,max(100,count(m)),km(m),(1.0E-9*cen(i,m),i=1,3) ! tcon 
!##################################################
!    Fragment 1 
!##################################################
   i=0; j=0; k=0 
   call ncell(tmp,i,j,k)
   if (debug(2) .eqv. .false.) then 
		write(mnum,'(i0,a,i0,a,i0)') m,"m",tmp,"r",ri ! 
		write(flname,'(a,a)') trim(mnum),'.com' 
		write(181,'(a,a)')  'g16  ',trim(flname) !, ".log"  !181 alone  
		call writegjf(flname,nam(m),c(:,:,m),elemt(:,m),ncpu,sets)
		call getmo(nam(m),elemt(:,m),mHOMO)
		write(mnum,'(2(a))') trim(mnum),".log" ! 
		if (levelm(m)) then 
			write(181,'(a)')  'orbital="error"'
			write(181,'(a,xa,xa)')  'orbital=`g09log', trim(mnum),'1  | tail -1`'
			write(181,'(a,xi0,xa)')  'echo', m , '$orbital  >> level.txt' 
			levelm(m)=.false.
		end if 
	end if 
   rt=rt+1   ! used to split run file

!##################################################
!    Fragment 2: neighbours 
!##################################################
   do l=1,max(100,count(m))
       read(131,*,end=350) tm,(frc(i),i=1,3),n !!!! cell 1 2 3 and molecules 
       if (tm /= m) goto 340 
       call ncell(tmp2,i,j,k); 

       ! if the dimer has already calculated then, not need to calculate again
       if (stt(m,n)) then
            write(flnm,'(2(i0,a))') m,'t',n,'.log' 
       else if(stt(n,m)) then 
            write(flnm,'(2(i0,a))') n,'t',m,'.log' 
       else
           !fragment 2 
           write(mnum2,'(i0,a,i0,a,i0)') n,"m",tmp2,"r",ri! 
           write(flname,'(a,a)') trim(mnum2),'.com' 

           if (debug(2) .eqv. .false.) write(181,'(a,a)')  'g16  ',trim(flname)  !,".log"  !181 alone 
             do iin=1,nam(n)
                dtc(1:3,iin)=c(1:3,iin,n)
                do i=1,3
                    dtc(1:3,iin)=dtc(1:3,iin)+frc(i)*cell(1:3,i)
                end do !i
            end do    
            
           if (debug(2) .eqv. .false.) call writegjf(flname,nam(n),dtc,elemt(:,n),ncpu,sets)
           call getmo(nam(n),elemt(:,n),nHOMO)
           write(mnum2,'(a,a)') trim(mnum2), ".log"! 
           if (levelm(n)) then ! initial is true
               if (debug(2) .eqv. .false.) write(181,'(a)')  'orbital="error"'
               if (debug(2) .eqv. .false.) write(181,'(a,xa,xa)')  'orbital=`g09log', trim(mnum2) ,'1  | tail -1`'
               if (debug(2) .eqv. .false.) write(181,'(a,xi0,xa)')  'echo', n , '$orbital  >> level.txt' 
               levelm(n)=.false.
           end if 
           rt=rt+1

!##################################################
!    Dimers
!##################################################
           write(flname,'(2(i0,a))') m,'t',n,'.com ' 
           if (debug(2) .eqv. .false.) write(181,'(a,a)')  'g16  ',flname !, ".log"  !181 dimer  
           write(flnm,'(2(i0,a))') m,'t',n,'.log' 
           if (debug(2) .eqv. .false.) write(181,'(a)')  'orbital="error"'
           if (debug(2) .eqv. .false.) write(181,'(a,xa,xa)')  'orbital=`g09log', trim(flnm) ,'1  | tail -1`'
           if (debug(2) .eqv. .false.) write(181,'(a,xa,xa)')  'echo', trim(flnm) , '$orbital  >> level.txt' 

           dc=0; delemt="NO"
           dc(:,1:nam(m))=c(:,1:nam(m),m);           
           delemt(1:nam(m))=elemt(1:nam(m),m)
           tchg(1:nam(m))=chg(1:nam(m))

           tmp=nam(m)+nam(n)
           dc(:,(nam(m)+1):tmp)=dtc(:,1:nam(n))
           delemt((nam(m)+1):tmp)=elemt(1:nam(n),n)
           tchg((nam(m)+1):tmp)=chg(1:nam(n))
   
           if (debug(2) .eqv. .false.) call writegjf(flname,tmp,dc,delemt,ncpu,sets) 
           if (debug(2) .eqv. .false.) write(181,"(a,xa20,xa20,xa20,xi0,xi0,xa,a,i0)")  & 
                    & "calv", mnum,mnum2,flnm, int(mhomo/2),int(nhomo/2),">>", "v.out",ri 
           if (debug(2) .eqv. .false.) write(181,"(3(xa))") "sleep 5; rm",mnum2,flnm
           stt(m,n)=.true.
           rt=rt+4
        end if 

        write(191,"(3(a))") "grep ' ",trim(flnm), " ' tmp.out >> v.out" 
   end do!n
340   backspace(131)
if (debug(2) .eqv. .false.) write(181,"(2(xa))") " rm ",mnum 
 
330 if (rt>32) then 
      close(181)
      ri=ri+1 
   end if 
end do !m

350 write(*,*) "end outputs"
close(131)
if (rt<=32) close(181)
close(191);close(211)
900  format(' ',i5,i5,*(xf9.5))
stop
end program ! main



subroutine ncell(n,i,j,k)
  implicit none
  integer::i,j,k,n
    n=(i+1)*3**2+(j+1)*3**1+(k+1)*3**0+1
  return 
end subroutine 
