program main 
!ccc molecule-atom-distance-basis-core cccccccccccccccccc
use lib 
implicit none 
integer::nam(10000),nm,i,j,k,l,m,n,count(10000),ii,rt,ri,tm,iin,iim  !number of atoms in molecule, and numbers of molecules
real::c(3,10000,10000),dc(3,10000),cen(3,10000),dis(3),cutoff,dtc(3,10000) !c(im,ia,3)
real::tdis,cell(3,3),kct,tdis1,tdis2,vtmp(3),dis1(3),laxs(3),tv1(3),tdis3
character(2)::elemt(10000,10000),delemt(20000) ! ,ALLOCATABLE
character(2),allocatable::tpielemt(:,:,:),pielemt(:,:,:)
character(20)::flname,mnum,mnum2,flnm,itpfile
integer::km(10000),ncpu,tmp,tmp2,piindx(4,100),jj,mindx(4),piindxnew(4,100) !,pair(1000,1000,2)
character(20)::sets 
character(200)::cmd
integer::nhomo,mhomo,rc(3),frc(3),nchain,chain(1000,10),ichain(10),tnv
logical::levelm(10000),stt(10000,10000)
real::morU(32),morM(32),morV(32),morP(32),morE(32),morC(32) ! 10000 molecule, one m has 10 chain, each change has 1000 atoms
real::tmp_pi(3,1000),sc
real,allocatable::picrd(:,:,:,:),tpicrd(:,:,:,:),cip(:,:),tpichg(:,:,:)
real::outparam(10),chg(1000),outparam1(10),outparam2(10),tchg(1000)
integer::kkm,kkn,nnc,nnca(100)
character(15)::method 
logical::debug(5),oldcon

debug=.False.

!debug(2)=.True.
method="snsr"

if (iargc()>=1) then
	CALL getarg(1, method)
    debug(2)=.True.
end if 

if (iargc()>=2) then 
	CaLL getarg(2, flnm)
	read(flnm,*) sc
else 
	sc=1.0
end if
!      write(*,*) 'Please input the coordinate file name'
!      read(*,*) flname
open(100,file="in",status='old') 
read(100,*) flname,itpfile
read(100,*) ncpu 
read(100,*) cutoff 
cutoff=cutoff*cutoff
read(100,*) 
read(100,"(a20)") sets
close(100)
!!!!!!!!!!!!!!!!!!!!!!!
!read itp file
!!!!!!!!!!!!


call readitp(itpfile,nchain,chain,ichain,chg,nnc,nnca)

if (debug(1)) then 
    i=1
    write(*,*) "number of atoms of chain 1: ", ichain(1)
    write(*,*) "chain: ",i
    do j=1,ichain(i)
        write(*,"(xi0)",ADVANCE="NO") chain(j,i)
    end do 
    write(*,*) 
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read molecules coordinate           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dc=0
call readgro(flname,cell,c,elemt,nm,nam,km) ! coordinate in nm
write(*,"(a,9(xf12.3))") "cell", cell
!Write(*,*) "Take the molecule separated by the boundary condition together" 
!write(*,*) "number of molecules",nm
do m=1,nm 
   do n=1,nam(m) 
     if (n/=1) then 
           rc=0
           do i=1,3
            if (c(i,n,m)-c(i,n-1,m)> 0.5*cell(i,i)) rc(i)=-1
            if (c(i,n,m)-c(i,n-1,m)<-0.5*cell(i,i)) rc(i)= 1
           end do !i
           do i=1,3
               c(1:3,n,m)=c(1:3,n,m)+rc(i)*cell(1:3,i)
           end do !i
     end if 
  end do ! n
end do ! m 

!do m=1,nm 
!   ri=0;tmp=0
!   write(mnum,'(i0,a,i0,a,i0)') m,"m",tmp,"r",ri ! 
!   write(flname,'(a,a)') trim(mnum),'.com' 
!   write(*,*) flname 
!   call writegjf(flname,nam(m),c(:,:,m),elemt(:,m),ncpu,sets)
!end do
!stop 

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



!##################################################
!           index for whole molecules    !
!##################################################
mindx=0
tdis1=0.d0
tdis2=0.d0
tdis3=0.d0
!use molecule 1 to find parameters

do j =1,nam(1)  ! find the longest axis in molecule
    do jj =j+1,nam(1)
        dis1(:)=c(:,jj,1)-c(:,j,1)
        tdis=dot(dis1(:),dis1(:)) 
        if (tdis>tdis1) then  ! find largest
            tdis1=tdis 
            mindx(1)=j
            mindx(2)=jj 
        end if 
    end do 
end do 

laxs=c(:,mindx(2),1)-c(:,mindx(1),1)
tnv=0

do i =1,nam(1)
    if (i .ne. mindx(1) .and. (i .ne. mindx(2))) then 
        if (tnv==0) then   !# whether we have set reference 
            tnv=1    !# set reference  
            vtmp=c(:,i,1)-c(:,mindx(1),1) ! #set as reference vector to distinguish left and right
        end if
        tv1=c(:,i,1)-c(:,mindx(1),1) !
        !theta=acos( dot(tv1,laxs) / sqrt(dot(laxs,laxs)*dot(tv1,tv1)) )    !#! cos= a*b/|a|b|
        !tdis=abs(dot(tv1,tv1)*sin(theta))
        tdis=sqrt(dot(cross(tv1,laxs),cross(tv1,laxs)))/sqrt(dot(laxs,laxs))  !|AXB|=|A||B|sin
        
        if ( dot( cross(vtmp,laxs),cross(tv1,laxs)) > 0.d00) then   !# left side 
            if (tdis>tdis2) then 
                tdis2=tdis !# find the largest dis 
                mindx(3)=i 
            end if
        else if (dot( cross(vtmp,laxs),cross(tv1,laxs)) < 0.d00 ) then !# right side 
            if (tdis>tdis3) then 
                tdis3=tdis !# find the largest dis 
                mindx(4)=i 
            end if
        end if
    end if
end do

write(*,"(a,4i5)") "molecule index", mindx
write(*,"(4(xa3,3(xf12.5)))")  elemt(mindx(1),1),c(:,mindx(1),1),elemt(mindx(2),1),c(:,mindx(2),1),&
& elemt(mindx(3),1),c(:,mindx(3),1),elemt(mindx(4),1),c(:,mindx(4),1)



!write(*,*) "pi index,npi", nchain
!if (debug(1)) write(*,*) "list all atoms of chain---1: ", chain(1:ichain(1),1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! use the first molecule to find 
!!!!! index atoms in pi planes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tmp=maxval(ichain(:))
allocate(picrd(3,tmp,nchain,nm),pielemt(tmp,nchain,nm),tpicrd(3,tmp,nchain,2),tpielemt(tmp,nchain,2),tpichg(tmp,nchain,2))

!picrd=0.d0
!do i=1,nm !   # find pi crd in all crd 
!    do k =1,nchain
!        do j =1,ichain(k)
!            picrd(:,j,k,i)=c(:,chain(j,k),i)
!        end do
!    end do
!end do

! index for pi planes
piindx=0
tnv=0
!#################
picrd=0.d0
do i=1,nm !   # find pi crd in all crd 
!i=1
    do k =1,nchain
        do j =1,ichain(k)
            picrd(:,j,k,i)=c(:,chain(j,k),i)
			pielemt(j,k,i)=elemt(chain(j,k),i)
        end do
    end do
end do

! index for pi planes
piindx=0
tnv=0

if (debug(1)) then 
    i=1
    write(*,*) "number of atoms of chain 1: ", ichain(1)
    write(*,*) "chain: ",i
    do j=1,ichain(i)
        write(*,"(xi0)",ADVANCE="NO") chain(j,i)
    end do 
    write(*,*) 
end if 

do k =1,nchain
    tmp_pi(:,1:ichain(k))=picrd(:,1:ichain(k),k,1)  ! (3,1000)# choose the first molecule to calculate index
	
   if (debug(2) .eqv. .false.) then 
		write(flname,'(a,i0,a)') '1c',k,'.com' 
		call writegjf(flname,ichain(k),picrd(:,1:ichain(k),k,1),pielemt(1:ichain(k),k,1),ncpu,sets)
	end if 
	
    tdis1=0.d0
    tdis2=0.d0
    tdis3=0.d0

    if (debug(1)) then 
        write(*,*) "chain: ", k, "atoms: ", ichain(k)
        write(*,*) "chain: ",k
        do j=1,ichain(i)
            write(*,"(xi0,3f8.4)") chain(j,k),tmp_pi(:,j)
        end do 
        write(*,*) 
    end if 

    do j =1,ichain(k)  ! find the longest atoms, three atoms of conjugate plane
        do jj =1,ichain(k)
            dis1(:)=tmp_pi(:,j)-tmp_pi(:,jj) 
            tdis=dot(dis1(:),dis1(:)) 
            if (tdis>tdis1) then
                tdis1=tdis 
                piindx(1,k)=chain(j,k)
                piindxnew(1,k)=j
                piindx(2,k)=chain(jj,k) 
                piindxnew(2,k)=jj
            end if 
        end do 
    end do 
    !  find atoms farest to the longest axis
    laxs=tmp_pi(:,piindxnew(1,k))-tmp_pi(:,piindxnew(2,k))
    
    !#########find the short axis
    tnv=0
    do i =1,ichain(k)
        if (i .ne. piindxnew(1,k) .and. i .ne. piindxnew(2,k)) then 
            if (tnv==0) then   !# whether we have set reference 
                tnv=1    !# set reference  
                vtmp=tmp_pi(:,i)-tmp_pi(:,piindxnew(1,k)) ! #set reference vector
                !vtmp=tmp_pi(:,i)-picrd(:,piindx(1,k),k,1)
                if (debug(1)) write(*,*) "set reference atom: ", chain(i,k) 
            else 
                tv1=tmp_pi(:,i)-tmp_pi(:,piindxnew(1,k)) !
                !tv1=tmp_pi(:,i)-picrd(:,piindx(1,k),k,1) 
                !theta=acos( dot(tv1,laxs) / sqrt(dot(laxs,laxs)*dot(tv1,tv1)) )    !#! cos= a*b/|a|b|
                !tdis=abs(dot(tv1,tv1)*sin(theta))
                tdis=sqrt(dot(cross(tv1,laxs),cross(tv1,laxs)))/sqrt(dot(laxs,laxs)) ! |AXB|=|A||B|sin, h=|A|sin=|AXB|/|B|
                
                if ( dot( cross(vtmp,laxs),cross(tv1,laxs)) > 0.d00) then   !# right side 
                    if (tdis>tdis2) then 
                        tdis2=tdis !# find the largest dis 
                        if (debug(1)) write(*,*) "right",chain(i,k),tdis2
                        piindx(3,k)=chain(i,k) 
                        piindxnew(3,k)=i 
                    end if
                else if (dot( cross(vtmp,laxs),cross(tv1,laxs)) < 0.d00 ) then !# left side 
                    if (tdis>tdis3) then 
                        tdis3=tdis !# find the largest dis 
                        if (debug(1)) write(*,*) "left",chain(i,k),tdis3
                        piindx(4,k)=chain(i,k) 
                        piindxnew(4,k)=i 
                    end if
                end if
            end if 
        end if
    end do
    write(*,"(a,4i5)") "k pi plane index", piindx(:,k)
    write(*,"(a,4i5)") "k pi plane new index", piindxnew(:,k)
    write(*,"(12(xf12.5))") (tmp_pi(:,piindxnew(i,k)),i=1,4)
end do 
!index finished 


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
    do n=1,nm
        tdis1=100
       if (n/=m ) then 
          do iim=1,nam(m)
          do iin=1,nam(n)
              vtmp(1:3)=c(1:3,iin,n)
              rc=0
              do i=1,3
               if (vtmp(i)-c(i,iim,m) >  0.5*cell(i,i)) rc(i)=-1
               if (vtmp(i)-c(i,iim,m) < -0.5*cell(i,i)) rc(i)= 1
              end do !i
              do i=1,3
                  vtmp(1:3)=vtmp(1:3)+rc(i)*cell(1:3,i)
              end do !i
              vtmp(1:3)=vtmp(1:3)-c(1:3,iim,m)
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
            if (tdis1 < 0.250) then 
                count(m)=count(m)+1
                vtmp(1:3)=cen(1:3,n)
                do i=1,3
                  vtmp(1:3)=vtmp(1:3)+frc(i)*cell(1:3,i)
                end do !i
                vtmp(1:3)=vtmp(1:3)-cen(1:3,m)  
               write(131,"(i5,3i3,i5,4(xf10.4))") m, frc(1:3),n,(vtmp(i),i=1,3),sqrt(tdis1) 
            end if 
         end if ! n/=m
    end do !n
!    write(*,*) "Finished", m 
!	if ((debug(2) .eqv. .true.) .and. (m >10) ) exit
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
write(cmd,"(a)") "m n S0 S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 &
           & S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25 S26 S27 S28 S29 S30 S31"
! output parameters
open(300,file="c3dmrs.dat",status='replace') 
open(301,file="u3dmrs.dat",status='replace') 
open(302,file="m3dmrs.dat",status='replace') 
open(303,file="v3dmrs.dat",status='replace') 
open(304,file="p3dmrs.dat",status='replace') 
open(305,file="e3dmrs.dat",status='replace') 
open(306,file="pigeomty.dat",status='replace') 
open(307,file="mgeomty.dat",status='replace') 
open(308,file="refpen.dat",status='replace') 
allocate(cip(nam(1),nam(1)))

open(400,file="dmrc3d.dat",status='replace') 
open(401,file="dmru3d.dat",status='replace') 
open(402,file="dmrm3d.dat",status='replace') 
open(403,file="dmrv3d.dat",status='replace') 
open(404,file="dmrp3d.dat",status='replace') 
open(405,file="dmre3d.dat",status='replace') 
open(406,file="rpi.dat",status='replace') 
open(407,file="rmol.dat",status='replace') 
open(408,file="dsall.dat",status='replace') 
open(409,file="neardisfft.dat",status='replace') 
open(410,file="cip.dat",status='replace') 
open(411,file="frame4.dat",status='replace') 


! output parameters
open(500,file="pic3d.dat",status='replace') 
open(501,file="piu3d.dat",status='replace') 
open(502,file="pim3d.dat",status='replace') 
open(503,file="piv3d.dat",status='replace') 
open(504,file="pip3d.dat",status='replace') 
open(505,file="pie3d.dat",status='replace') 

open(600,file="dpic3d.dat",status='replace') 
open(601,file="dpiu3d.dat",status='replace') 
open(602,file="dpim3d.dat",status='replace') 
open(603,file="dpiv3d.dat",status='replace') 
open(604,file="dpip3d.dat",status='replace') 
open(605,file="dpie3d.dat",status='replace') 
open(606,file="picip.dat",status='replace') 

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

write(*,*) "begin output g09 and descriptors"
do m=1,nm 
!##################################################
!    Generate list of job need to run 
!##################################################
  if (rt>32) then   ! used to split run file, count number of calculations
      write(flname,"(a,i0)") "run",ri 
      if (debug(2) .eqv. .false.) open(181, file =flname, status = 'replace') !command to run g09
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
		write(181,'(a,a)')  'g09  ',trim(flname) !, ".log"  !181 alone  
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

           if (debug(2) .eqv. .false.) write(181,'(a,a)')  'g09  ',trim(flname)  !,".log"  !181 alone 
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
           if (debug(2) .eqv. .false.) write(181,'(a,a)')  'g09  ',flname !, ".log"  !181 dimer  
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
!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate 3D MORSE code
!!!!!!!!!!!!!!!!!!!!!!!!!
            call morse3d(sc,method,dc(1:3,1:tmp),delemt(1:tmp),tmp,morU,morM,morV,morP,morE,tchg,morC) !normal 3D-MORSE
            write(300,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morC
            write(301,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morU
            write(302,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morM
            write(303,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morV
            write(304,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morP
            write(305,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morE

			call morse3ddimer(sc,method,c(:,1:nam(m),m),elemt(1:nam(m),m),nam(m),&      !CIP,  exclude intramolecular interaction
					&dtc(:,1:nam(n)),elemt(1:nam(n),n),nam(n),&
					& morU,morM,morV,morP,morE,tchg,tchg,morC,cip)
            write(400,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morC
            write(401,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morU
            write(402,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morM
            write(403,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morV
            write(404,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morP
            write(405,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morE
			write(410,900) m,n,cip(1:nam(m),1:nam(n))

! distance of 4v, output 16 distance 
            call  calparam4(mindx(1:4),dc(:,1:nam(m)),nam(m),mindx(1:4),dtc(:,1:nam(n)),nam(n),morU) 
            write(411,"(i5,i5,16(xf13.5))") m,n,(morU(i),i=1,16)

			call dnearest(c(:,1:nam(m),m),elemt(1:nam(m),m),nam(m),&    ! !  nearest 32th distance 
					&dtc(:,1:nam(n)),elemt(1:nam(n),n),nam(n),&
					& morU,morM,morV,morP,morE,tchg,tchg,morC)
			write(408,"(i5,i5,12(xf13.5),20(xf11.5))",ADVANCE='NO') m,n, (morU(i),i=1,32)
            call calparam4v(mindx(1:4),dc(:,1:nam(m)),nam(m),&        !! OUTFRAME OF MOLECULES
            mindx(1:4),dtc(:,1:nam(n)),nam(n),morU)
            write(408,"(22(xf10.5))") (morU(i),i=1,22) 

			call msnearest(c(:,1:nam(m),m),elemt(1:nam(m),m),nam(m),&    !! 3D-MORSE calculated on nearest atoms
					&dtc(:,1:nam(n)),elemt(1:nam(n),n),nam(n),&
					& morU,morM,morV,morP,morE,tchg,tchg,morC)
			write(409,"(i5,i5,12(xf13.5),20(xf11.5))") m,n,morU ! 32

!!!!!descriptors: seperated whole molecule                              
            call calparam2(mindx(1:4),dc(:,1:nam(m)),nam(m),&           !!!!!!!!! MOL 
            mindx(1:4),dtc(:,1:nam(n)),nam(n),morU)
            write(307,"(2i5,20(xf10.5),xf10.5)") m,n,(morU(i),i=1,21)

            call calparam2r(mindx(1:4),dc(:,1:nam(m)),nam(m),&     !!!!!!!!! MOL, with accurate height of pi plane,but not good
            mindx(1:4),dtc(:,1:nam(n)),nam(n),morU)
            write(407,"(2i5,20(xf10.5),xf10.5)") m,n,(morU(i),i=1,21)

            call calref(mindx(1:4),nam(m),dc(:,1:nam(m)),delemt(1:nam(m)),&    !!!Lederer's
            mindx(1:4),nam(n),dtc(:,1:nam(n)),delemt(1:nam(m)),outparam)
            write(308,"(2i5,6(xf10.5))") m,n,outparam(1:6)

!!!!!! write descriptor of pi planes (independent of position)!!!!!!!!!!!
            do kkm =1,nchain
                do ii =1,ichain(kkm)
                    !picrd(:,ii,kkm,m)=c(:,chain(ii,kkm),m)
                    tpicrd(:,ii,kkm,1)=dc(:,chain(ii,kkm))
                    tpielemt(ii,kkm,1)=elemt(chain(ii,kkm),m) 
                    tpichg(ii,kkm,1)=chg(chain(ii,kkm))  
                end do 
            end do
            do kkn =1,nchain
                do ii =1,ichain(kkn)
                    !picrd(:,ii,kkn,n)=c(:,chain(ii,kkn),n)
                    tpicrd(:,ii,kkn,2)=dtc(:,chain(ii,kkn))
                    tpielemt(ii,kkn,2)=elemt(chain(ii,kkn),n) 
                    tpichg(ii,kkn,2)=chg(chain(ii,kkn)) 
                end do 
            end do

!!!!!descriptors: seperated pi planes 
            write(306,"(2i5)",ADVANCE='NO') m,n
            do kkm =1,nchain
                call calparam1(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),outparam)   !! mono molecular
                write(306,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
                call calparam1r(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),outparam)  !! revised mono molecular
                write(406,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
            end do 
            do kkn =1,nchain
                call calparam1(piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),outparam)   !! mono molecular
                write(306,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
                call calparam1r(piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),outparam)  !! revised mono molecular
                write(406,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
            end do 
			!write(*,*) "stage 1";read(*,*) kkn
!!!!! descriptors:  pi plane-- pi plane
                   write(500,"(i5,i5)",advance='no') m,n
                   write(501,"(i5,i5)",advance='no') m,n
                   write(502,"(i5,i5)",advance='no') m,n
                   write(503,"(i5,i5)",advance='no') m,n
                   write(504,"(i5,i5)",advance='no') m,n
                   write(505,"(i5,i5)",advance='no') m,n
                   write(600,"(i5,i5)",advance='no') m,n
                   write(601,"(i5,i5)",advance='no') m,n
                   write(602,"(i5,i5)",advance='no') m,n
                   write(603,"(i5,i5)",advance='no') m,n
                   write(604,"(i5,i5)",advance='no') m,n
                   write(605,"(i5,i5)",advance='no') m,n
                   write(606,900,ADVANCE='NO') m,n
			!write(*,*) "stage 2";read(*,*) kkn
            do kkm=1,nchain
                do kkn=1,nchain 
                    call calparam2(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),&    !  PI
                    & piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),morU)
                    !write(306,"(8(xf10.5),xf10.5)",ADVANCE='NO') outparam(1:9) 
                    write(306,"(9(xf10.5))",ADVANCE='NO') (morU(i),i=13,21)
                    call calparam2r(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),&   ! REVISED PI
                    & piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),morU)
                    !write(306,"(8(xf10.5),xf10.5)",ADVANCE='NO') outparam(1:9) 
                    write(406,"(9(xf10.5))",ADVANCE='NO') (morU(i),i=13,21)
				!write(*,*) "stage 2.1";read(*,*) i
                   dc=0; delemt="NO"
                   dc(:,1:ichain(kkm))=tpicrd(:,:,kkm,1)           
                   delemt(1:ichain(kkm))=tpielemt(:,kkm,1)
                   tchg(1:ichain(kkm))=tpichg(ii,kkm,1)

                   tmp=ichain(kkm)+ichain(kkn)
                   dc(:,(ichain(kkm)+1):tmp)=tpicrd(:,:,kkn,2)
                   delemt((ichain(kkm)+1):tmp)=tpielemt(:,kkn,2)
                   tchg((ichain(kkm)+1):tmp)=chg(1:nam(n))
   				!write(*,*) "stage 2.2";read(*,*) i
                   call morse3d(sc,method,dc(1:3,1:tmp),delemt(1:tmp),tmp,morU,morM,morV,morP,morE,tchg,morC) !normal 3D-MORSE
                   write(500,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morC
                   write(501,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morU
                   write(502,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morM
                   write(503,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morV
                   write(504,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morP
                   write(505,"(12(xf13.5),20(xf11.5))",ADVANCE='NO') morE
				!write(*,*) "stage 2.3";read(*,*) i
				
                   call morse3ddimer(sc,method,tpicrd(:,:,kkm,1),tpielemt(:,kkm,1),ichain(kkm),&      !CIP,  exclude intramolecular interaction
                           &tpicrd(:,:,kkn,2),tpielemt(:,kkn,2),ichain(kkn),&
                           & morU,morM,morV,morP,morE,tchg,tchg,morC,cip)
                   write(600,"(12(xf13.5),20(xf11.5))",advance='no') morC
                   write(601,"(12(xf13.5),20(xf11.5))",advance='no') morU
                   write(602,"(12(xf13.5),20(xf11.5))",advance='no') morM
                   write(603,"(12(xf13.5),20(xf11.5))",advance='no') morV
                   write(604,"(12(xf13.5),20(xf11.5))",advance='no') morP
                   write(605,"(12(xf13.5),20(xf11.5))",advance='no') morE
				!write(*,*) "stage 2.4";read(*,*) i
                   write(606,'(*(xf9.5))',advance='no') cip(1:ichain(kkm),1:ichain(kkn))
                end do 
            end do
			!write(*,*) "stage 3"
			!read(*,*) kkn
            write(306,*);write(406,*)
			write(500,*);write(501,*);write(502,*);write(503,*);write(504,*);write(505,*)
			write(600,*);write(601,*);write(602,*);write(603,*);write(604,*);write(605,*);write(606,*)
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
close(400);close(401);close(402);close(403);close(404);close(405);close(406);close(407);close(408);close(409)
close(410);close(411)
close(300);close(301);close(302);close(303);close(304);close(305);close(306);close(307);close(308)
close(600);close(601);close(602);close(603);close(604);close(605);close(606)
close(500);close(501);close(502);close(503);close(504);close(505)
stop 

!340 write(*,*) "tm is not equal to m", tm,m 
close(131)
if (rt<=32) close(181)
close(191)
close(211) 
stop
900  format(' ',i5,i5,*(xf9.5))
end program ! main



subroutine ncell(n,i,j,k)
  implicit none
  integer::i,j,k,n
    n=(i+1)*3**2+(j+1)*3**1+(k+1)*3**0+1
  return 
end subroutine 
