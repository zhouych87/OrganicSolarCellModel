program main 
!ccc molecule-atom-distance-basis-core cccccccccccccccccc
use lib 
use MPI 
implicit none 
integer::nam(10000),nm,i,j,k,l,m,n,count(1000),ii,rt,ri,tm,iin,iim  !number of atoms in molecule, and numbers of molecules
real::c(3,1000,2),dc(3,1000),cen(3,1000),dis(3),cutoff,dtc(3,1000) !c(im,ia,3)
real::tdis,cell(3,3),kct,tdis1,tdis2,vtmp(3),dis1(3),laxs(3),tv1(3),tdis3
character(2)::elemt(1000,2),delemt(2000) ! ,ALLOCATABLE
character(2),allocatable::tpielemt(:,:,:),pielemt(:,:,:)
character(20)::flname,mnum,mnum2,flnm,itpfile
integer::km(10000),ncpu,tmp,tmp2,piindx(4,100),jj,mindx(4),piindxnew(4,100) !,pair(1000,1000,2)
character(20)::sets 
character(6)::label(1000)
character(200)::cmd
integer::nhomo,mhomo,rc(3),frc(3),frc0(3),nchain,chain(1000,10),ichain(10),tnv,rc0(3)
real,allocatable::morU(:),morM(:),morV(:),morP(:),morE(:),morC(:) ! 10000 molecule, one m has 10 chain, each change has 1000 atoms
real::tmp_pi(3,1000),sc
real,allocatable::picrd(:,:,:,:),tpicrd(:,:,:,:),cip(:,:),tpichg(:,:,:)
real::outparam(10),chg(1000),outparam1(10),outparam2(10),tchg(1000),tchg1(1000)
integer::anglhist(100),rhist(100),rchist(100)
integer::kkm,kkn,nnc,nnca(100),an,stat,nsc,m0,n0
character(15)::method,mform
logical::debug(5),oldcon,tindex,grof,voutf
character(20)::grofn
real::allc(3,1000,10000)
character(2)::allelemt(1000,10000)
integer :: myid,ierr,npcs,status(MPI_STATUS_SIZE)

call MPI_INIT( ierr )     
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, npcs, ierr )

debug=.False.
tindex=.false.
grof = .false.

if (myid .eq. 0) then 
    if (iargc()==0) then
        write(*,*) "using  : desgen itpfile method mform coefficient grofn"
        write(*,*) "methods: geo, 3dm, d3d, 3dpi, d3dpi"
        write(*,*) "mform  : snsr,cnsr,snr, cnr, snsqrtr,cnsqrtr,sn,cn, and their abs: abssnsr ...."
        write(*,*) "coefficient in 3dm"
        stop
    end if
end if 

CALL getarg(1, itpfile)
CALL getarg(2, method)
if (iargc()>=3) then 
    call getarg(3, mform)
else 
    mform="snsr"
end if 
if (iargc()>=4) then
    call getarg(4, flnm)
    read(flnm,*) sc
else 
    sc=1.0 
end if
if (iargc()>=5) call getarg(5,grofn)
if (iargc()==6) then 
    call getarg(6,flnm)
    read(flnm,*) nsc
else 
    nsc=32
end if 

allocate(morU(nsc),morM(nsc),morV(nsc),morP(nsc),morE(nsc),morC(nsc))

call readitp(itpfile,nchain,chain,ichain,chg,nnc,nnca,an)
if (debug(1)) then
    i=1
    if (myid .eq. 0) write(*,*) "number of atoms of chain 1: ", ichain(1)
    if (myid .eq. 0) write(*,*) "chain: ",i
    do j=1,ichain(i)
        if (myid .eq. 0) write(*,"(xi0)",ADVANCE="NO") chain(j,i)
    end do
    if (myid .eq. 0) write(*,*)
end if


inquire(file="v.out", exist=voutf)
if (voutf .eqv. .true.) then 
    if (myid .eq. 0) print *, "Will use v.out to read dimer information"
    open(101,file="v.out",status='old')
else 
    inquire(file="connection.dat", exist=oldcon)
    if (oldcon .eqv. .false.) then
        if (myid .eq. 0) write(*,*) "no connection.dat is found. Will generate com files"
    else !!!!! connection.dat exist
        if (myid .eq. 0) write(*,*) "connection.dat exist, no com and run file will be writen"
        open(102,file="connection.dat",status='old')
    end if 
end if

if (iargc()>=5) then !caLL getarg(5,grofn)
    inquire(file=grofn, exist=grof)
    if (grof .eqv. .false.) then
        if (myid .eq. 0) print *, grofn, " gro file is not exist" 
        !stop
    else 
        call readgro(grofn,cell,allc,allelemt,nm,nam,km)
        cell=10.0*cell;allc=10.0*allc ! from nm to A
        if (myid .eq. 0) write(*,"(a,9(xf12.3))") "cell", cell
        do m=1,nm 
           do n=1,nam(m) 
             if (n/=1) then 
                   rc=0
                   do i=1,3
                    if (allc(i,n,m)-allc(i,n-1,m)> 0.5*cell(i,i)) rc(i)=-1
                    if (allc(i,n,m)-allc(i,n-1,m)<-0.5*cell(i,i)) rc(i)= 1
                   end do !i
                   do i=1,3
                       allc(1:3,n,m)=allc(1:3,n,m)+rc(i)*cell(1:3,i)
                   end do !i
             end if 
          end do ! n
        end do ! m 
        print *, "read gro completed"
    end if 
end if


!##################################################
!    Generate descriptor files
!##################################################
!write(cmd,"(a)") "m n S0 S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 &
!           & S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25 S26 S27 S28 S29 S30 S31"
! output parameters
allocate(cip(an,an))
if (myid .eq. 0) then 
    if (trim(method) == 'geo') then
        write(*,*) "Generate geometry based descriptors: open file "
        open(307,file="zhougeo.dat",status='replace')
        open(308,file="geo.dat",status='replace')
        open(309,file="angle3d.dat",status='replace')
        open(310,file="rdfhist.dat",status='replace',iostat=stat)
        open(407,file="rmol.dat",status='replace')
        open(408,file="dsall.dat",status='replace')
        open(409,file="neardisfft.dat",status='replace')
        open(410,file="cip.dat",status='replace')
        open(411,file="frame4.dat",status='replace')
    else if (trim(method) == 'geopi') then
        write(*,*) "Generate geometry based descriptors: open file "
        open(306,file="pizhougeopi.dat",status='replace')
        open(406,file="pir.dat",status='replace')
        open(606,file="picip.dat",status='replace')
        open(609,file="piangle3d.dat",status='replace')
        open(610,file="pirdfhist.dat",status='replace')
    else if (trim(method) == '3dm') then
        write(*,*) "Generate 3D MORSE descriptors"
        open(300,file="c3dm.dat",status='replace')
        open(301,file="u3dm.dat",status='replace')
        open(302,file="m3dm.dat",status='replace')
        open(303,file="v3dm.dat",status='replace')
        open(304,file="p3dm.dat",status='replace')
        open(305,file="e3dm.dat",status='replace')
    else if (trim(method) == 'd3d') then
        write(*,*) "Generate intermolecular 3D MORSE descriptors"
        open(400,file="dc3d.dat",status='replace')
        open(401,file="du3d.dat",status='replace')
        open(402,file="dm3d.dat",status='replace')
        open(403,file="dv3d.dat",status='replace')
        open(404,file="dp3d.dat",status='replace')
        open(405,file="de3d.dat",status='replace')
    else if (trim(method) == '3dpi') then
        write(*,*) "Generate pi 3D MORSE descriptors"
        open(500,file="c3dpi.dat",status='replace')
        open(501,file="u3dpi.dat",status='replace')
        open(502,file="m3dpi.dat",status='replace')
        open(503,file="v3dpi.dat",status='replace')
        open(504,file="p3dpi.dat",status='replace')
        open(505,file="e3dpi.dat",status='replace')
    else if (trim(method) == 'd3dpi') then
        write(*,*) "Generate intermolecular pi 3D MORSE descriptors"
        open(600,file="dc3dpi.dat",status='replace')
        open(601,file="du3dpi.dat",status='replace')
        open(602,file="dm3dpi.dat",status='replace')
        open(603,file="dv3dpi.dat",status='replace')
        open(604,file="dp3dpi.dat",status='replace')
        open(605,file="de3dpi.dat",status='replace')
    end if 
end if
!!!!!!!!!!!!!!!!!!!!!!
!!!!LOOP!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
do while (.true.)
    if (voutf .eqv. .true.)  then 
        !write(*,*) "read v.out"
        read(101,*,end=910) m,n ! v.out
    else 
        if (oldcon .eqv. .true.) then
            read(102,*,end=910) m,(frc0(i),i=1,3),n !connection.dat
        end if 
    end if 

    if (grof .eqv. .false.) then
    !if (iargc()<5) then     ! no gro is provided, read com or gjf file to generate descriptors
        write(flname,"(i0,a,i0,a)") m,"t",n,".gjf"; inquire(file=flname, exist=oldcon)
        if (oldcon .eqv. .false.) then
            write(flname,"(i0,a,i0,a)") n,"t",m,".gjf"; inquire(file=flname, exist=oldcon)
            if (oldcon .eqv. .false.) then
                write(flname,"(i0,a,i0,a)") n,"t",m,".com"; inquire(file=flname, exist=oldcon)
                if (oldcon .eqv. .false.) then
                    write(flname,"(i0,a,i0,a)") m,"t",n,".com"; inquire(file=flname, exist=oldcon)
                    if (oldcon .eqv. .false.) then
                        if (myid .eq. 0) write(*,*) flname, "and .com is not found"
                        stop 
                    end if
                end if
            end if 
        end if 
        if (myid .eq. 0) write(*,*) 'Will read ',flname
        call readgjf(flname,an,c,elemt)
    else 
        !print *, "an: ",nam(1),m,n 
        an=nam(m)
        c(:,1:an,1)=1.d0*allc(:,1:an,m)  ! 
        c(:,1:an,2)=1.d0*allc(:,1:an,n)

        tdis1=100.0
        do iim=1,an
        do iin=1,an
            vtmp(1:3)=c(1:3,iin,1)
            rc=0
            do i=1,3
             if (vtmp(i)-c(i,iim,2) >  0.5*cell(i,i)) rc(i)= 1 ! vtmp is too big
             if (vtmp(i)-c(i,iim,2) < -0.5*cell(i,i)) rc(i)=-1
            end do !i
            do i=1,3
                vtmp(1:3)=vtmp(1:3)-rc(i)*cell(1:3,i) ! has to be deducted.
            end do !i
            vtmp(1:3)=vtmp(1:3)-c(1:3,iim,2)

            tdis=vtmp(1)*vtmp(1)+vtmp(2)*vtmp(2)+vtmp(3)*vtmp(3)
            if (tdis<tdis1) then 
              tdis1=tdis  ! find minimum to tdis1 
              frc=rc
            end if 
        end do !iin
        end do !iim

!        if (frc0(1) .ne. frc(1) .or. frc0(2) .ne. frc(2) .or. frc0(3) .ne. frc(3) &
!   & .or. m.ne.m0 .or. n .ne. n0 ) then 
!           write(*,*) m,n, frc, m0,n0, frc0
!           stop 
 !       end if
        do iin=1,an  ! for across boundary dimers , n is the second
            dtc(1:3,iin)=1.d0*allc(:,iin,n) 
            do i=1,3
                dtc(1:3,iin)=dtc(1:3,iin)+frc(i)*cell(1:3,i)*1.d0
            end do !i
        end do 

        c(:,1:an,2)=dtc(1:3,1:an)
        elemt(1:an,1)=allelemt(1:an,m)
        elemt(1:an,2)=allelemt(1:an,n)
    end if 

    !write(*,*) "start first molecule"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! first dimer will be used for index
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (tindex .eqv. .false.) then 
        if (myid .eq. 0) write(*,*) "The first dimer are molecule",n,m
        if (myid .eq. 0) write(*,*) "Start indexing"
        tindex=.true.
        !geo
        do l=1,2
            do k=1,an
                cen(:,l)=cen(:,l)+c(:,k,l)
            end do
            cen(:,l)=cen(:,l)/an
        end do
        
        mindx=0
        tdis1=0.d0
        tdis2=0.d0
        tdis3=0.d0
        !use molecule 1 to find parameters
        do j =1,an  ! find the longest axis in molecule
            do jj =j+1,an
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
        do i =1,an
            if (i .ne. mindx(1) .and. (i .ne. mindx(2))) then
                if (tnv==0) then   !# whether we have set reference
                    tnv=1    !# set reference
                    vtmp=c(:,i,1)-c(:,mindx(1),1) ! #set as reference vector to distinguish left and right
                end if
                tv1=c(:,i,1)-c(:,mindx(1),1) !
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
        
        if (myid .eq. 0) write(*,"(a,4i5)") "molecule index", mindx
        if (myid .eq. 0) write(*,"(4(xa3,3(xf12.5)))")  elemt(mindx(1),1),c(:,mindx(1),1),elemt(mindx(2),1),c(:,mindx(2),1),&
        & elemt(mindx(3),1),c(:,mindx(3),1),elemt(mindx(4),1),c(:,mindx(4),1)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! find index atoms in pi planes
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        tmp=maxval(ichain(:))
        allocate(picrd(3,tmp,nchain,2),pielemt(tmp,nchain,2),tpicrd(3,tmp,nchain,2),tpielemt(tmp,nchain,2),tpichg(tmp,nchain,2))
        ! index for pi planes
        piindx=0
        tnv=0
        !#################
        picrd=0.d0
        do i=1,2 !   # find pi crd in all crd
        !i=1
            do k =1,nchain
                do j =1,ichain(k)
                    picrd(:,j,k,i)=c(:,chain(j,k),i)   ! has to be label
                    pielemt(j,k,i)=elemt(chain(j,k),i)
                end do
            end do
        end do
    end if 
!!!!!!!!!!!!!!!!!!!!!!!
!  end index !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!
    !if (m==1) write(*,*) "begin crd given"
    dc=0; delemt="NO"
    dc(  :,1:an)=c(:  ,1:an,1);
    delemt(1:an)=elemt(1:an,1)
    tchg(  1:an)=  chg(1:an)
    
    tmp=an+an
    dc(  :,(an+1):tmp)=  c(:,1:an,2);
    delemt((an+1):tmp)=elemt(1:an,2)
    tchg(  (an+1):tmp)=  chg(1:an)
    
    if (myid .eq. 0 .and. m==1) write(*,*) "begin output descriptors"
    if (trim(method) == '3dm') then
        call morse3d(sc,mform,dc(1:3,1:tmp),delemt(1:tmp),tmp,morU,morM,morV,morP,morE,morC,tchg,nsc,npcs,myid,ierr) !normal 3D-MORSE
        if (myid .eq. 0) write(300,901) m,n,morC
        if (myid .eq. 0) write(301,901) m,n,morU
        if (myid .eq. 0) write(302,901) m,n,morM
        if (myid .eq. 0) write(303,901) m,n,morV
        if (myid .eq. 0) write(304,901) m,n,morP
        if (myid .eq. 0) write(305,901) m,n,morE
        !write(*,*) m,n, "1 heading"
    else if (trim(method) == 'd3d') then
        call morse3ddimer(sc,mform,c(:,1:an,1),elemt(1:an,1),an,&      !exclude intramolecular interaction
                &c(:,1:an,2),elemt(1:an,2),an,&
                & morU,morM,morV,morP,morE,morC,chg,chg,nsc,npcs,myid,ierr)
        if (myid .eq. 0) write(400,901) m,n,morC
        if (myid .eq. 0) write(401,901) m,n,morU
        if (myid .eq. 0) write(402,901) m,n,morM
        if (myid .eq. 0) write(403,901) m,n,morV
        if (myid .eq. 0) write(404,901) m,n,morP
        if (myid .eq. 0) write(405,901) m,n,morE
    
    else if (trim(method) == 'geo') then
        call calcip(an,c(:,1:an,1),elemt(1:an,1),&
            &an, c(:,1:an,2),elemt(1:an,2),chg(1:an),chg(1:an),cip(1:an,1:an),anglhist,rhist,npcs,myid,ierr)
        if (myid .eq. 0) write(410,900) m,n,cip(1:an,1:an)
        if (myid .eq. 0) write(310,"(*(i5))") m,n,rhist !310
        if (myid .eq. 0) write(309,"(*(i5))")  m,n,anglhist !309
    
    ! distance of 4v, output 16 distance
        call  calparam4(mindx(1:4),c(:,1:an,1),an,mindx(1:4),c(:,1:an,2),an,morU(1:16))
        if (myid .eq. 0) write(411,"(i5,i5,16(xf13.5))") m,n,(morU(i),i=1,16)
        call dnearest(c(:,1:an,1),elemt(1:an,1),an,&    ! !  nearest 32th distance
            &c(:,1:an,2),elemt(1:an,2),an,&
            & morU(1:32))
        if (myid .eq. 0) write(408,"(i5,i5,12(xf13.5),20(xf11.5))",ADVANCE='NO') m,n, (morU(i),i=1,32)
        call calparam4v(mindx(1:4),c(:,1:an,1),an,&        !! OUTFRAME OF MOLECULES
        mindx(1:4),c(:,1:an,2),an,morU)
        if (myid .eq. 0) write(408,"(22(xf10.5))") (morU(i),i=1,22)

    !!!!!descriptors: seperated whole molecule
        call calparam2(mindx(1:4),c(:,1:an,1),an,&           !!!!!!!!! MOL
        mindx(1:4),c(:,1:an,2),an,morU)
        if (myid .eq. 0) write(307,"(2i5,20(xf10.5),xf10.5)") m,n,(morU(i),i=1,21)
        call calparam2r(mindx(1:4),c(:,1:an,1),an,&     !!!!!!!!! MOL, with accurate height of pi plane,but not good
        mindx(1:4),c(:,1:an,2),an,morU)
        if (myid .eq. 0) write(407,"(2i5,20(xf10.5),xf10.5)") m,n,(morU(i),i=1,21)
        call calref(mindx(1:4),an,c(:,1:an,1),delemt(1:an),&    !!!Lederer's
        mindx(1:4),an,c(:,1:an,2),delemt(1:an),outparam)
        if (myid .eq. 0) write(308,"(2i5,6(xf10.5))") m,n,outparam(1:6)
    
    else if (trim(method) == 'd3dpi' .or. trim(method) == '3dpi'.or. trim(method) == 'geopi') then
    !!!!!! write descriptor of pi planes (independent of position)!!!!!!!!!!!
        do kkm =1,nchain
            do ii =1,ichain(kkm)
                tpicrd(:,ii,kkm,1)=   c(:,chain(ii,kkm),1)  !dc(:,chain(ii,kkm))
                tpielemt(ii,kkm,1)= elemt(chain(ii,kkm),1)
                tpichg(  ii,kkm,1)=   chg(chain(ii,kkm))
            end do
        end do
        do kkn =1,nchain
            do ii =1,ichain(kkn)
                tpicrd(:,ii,kkn,2)=   c(:,chain(ii,kkn),2)
                tpielemt(ii,kkn,2)= elemt(chain(ii,kkn),2)
                tpichg(  ii,kkn,2)  = chg(chain(ii,kkn))
            end do
        end do
    
    !!!geopi
        if (trim(method) == 'geopi') then
    !!!!!descriptors: seperated pi planes
            if (myid .eq. 0) write(306,"(2i5)",ADVANCE='NO') m,n
            do kkm =1,nchain
                call calparam1(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),outparam)   !! mono molecular
                if (myid .eq. 0) write(306,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
                call calparam1r(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),outparam)  !! revised mono molecular
                if (myid .eq. 0) write(406,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
            end do
            do kkn =1,nchain
                call calparam1(piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),outparam)   !! mono molecular
                if (myid .eq. 0) write(306,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
                call calparam1r(piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),outparam)  !! revised mono molecular
                if (myid .eq. 0) write(406,"(6(xf10.5))",ADVANCE='NO') outparam(1:6)
            end do
            if (myid .eq. 0) write(606,900,ADVANCE='NO') m,n
            !write(*,*) "stage 2";read(*,*) kkn
            do kkm=1,nchain
                do kkn=1,nchain
                    call calparam2(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),&    !  PI
                    & piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),morU)
                    !write(306,"(8(xf10.5),xf10.5)",ADVANCE='NO') outparam(1:9)
                    if (myid .eq. 0) write(306,"(9(xf10.5))",ADVANCE='NO') (morU(i),i=13,21)
                    call calparam2r(piindxnew(1:4,kkm),tpicrd(:,:,kkm,1),ichain(kkm),&   ! REVISED PI
                    & piindxnew(1:4,kkn),tpicrd(:,:,kkn,2),ichain(kkn),morU)
                    !write(306,"(8(xf10.5),xf10.5)",ADVANCE='NO') outparam(1:9)
                    if (myid .eq. 0) write(406,"(9(xf10.5))",ADVANCE='NO') (morU(i),i=13,21)
        
                    call calcip(ichain(kkm),tpicrd(:,:,kkm,1),tpielemt(:,kkm,1),&      !CIP,  exclude intramolecular interaction
                            &ichain(kkn),tpicrd(:,:,kkn,2),tpielemt(:,kkn,2),&
                            & tpichg(:,kkm,1),tpichg(:,kkn,2),cip,anglhist,rhist,npcs,myid,ierr)
                    if (myid .eq. 0) write(606,900) m,n,cip
                    if (myid .eq. 0) write(609,"(100(i5))",advance='no') anglhist
                    if (myid .eq. 0) write(610,"(100(i5))",advance='no') rhist
                    !write(*,*) m,n, "1 ending"
                end do
            end do
            if (myid .eq. 0) then 
                write(306,*);write(406,*);write(606,*);write(609,*);write(610,*)
            end if 
    
        else if (trim(method) == '3dpi') then
    !!!!! descriptors:  pi plane-- pi plane
            if (myid .eq. 0) then 
            write(500,"(i5,i5)",advance='no') m,n
            write(501,"(i5,i5)",advance='no') m,n
            write(502,"(i5,i5)",advance='no') m,n
            write(503,"(i5,i5)",advance='no') m,n
            write(504,"(i5,i5)",advance='no') m,n
            write(505,"(i5,i5)",advance='no') m,n
            end if 
            do kkm=1,nchain
                do kkn=1,nchain
            
                dc=0; delemt="NO"
                dc(:  ,1:ichain(kkm))=tpicrd(:,:,kkm,1)
                delemt(1:ichain(kkm))=tpielemt(:,kkm,1)
                tchg(  1:ichain(kkm))=  tpichg(:,kkm,1)
            
                tmp=ichain(kkm)+ichain(kkn)
                dc(:,  (ichain(kkm)+1):tmp)=tpicrd(:,:,kkn,2)
                delemt((ichain(kkm)+1):tmp)=tpielemt(:,kkn,2)
                tchg(  (ichain(kkm)+1):tmp)=  tpichg(:,kkn,2)
            
                call morse3d(sc,mform,dc(1:3,1:tmp),delemt(1:tmp),tmp,morU,morM,morV,morP,morE,morC,tchg(1:tmp),nsc,npcs,myid,ierr) !normal 3D-MORSE
                if (myid .eq. 0) then 
                write(500,902,ADVANCE='NO') morC
                write(501,902,ADVANCE='NO') morU
                write(502,902,ADVANCE='NO') morM
                write(503,902,ADVANCE='NO') morV
                write(504,902,ADVANCE='NO') morP
                write(505,902,ADVANCE='NO') morE
                end if 
            end do
            end do
            if (myid .eq. 0) then 
                write(500,*);write(501,*);write(502,*);write(503,*);write(504,*);write(505,*)
            end if 
        
        else if ( trim(method) == 'd3dpi') then 
        if (myid .eq. 0) then 
            write(600,"(i5,i5)",advance='no') m,n
            write(601,"(i5,i5)",advance='no') m,n
            write(602,"(i5,i5)",advance='no') m,n
            write(603,"(i5,i5)",advance='no') m,n
            write(604,"(i5,i5)",advance='no') m,n
            write(605,"(i5,i5)",advance='no') m,n
        end if 
            do kkm=1,nchain
                do kkn=1,nchain
                    call morse3ddimer(sc,mform,tpicrd(:,:,kkm,1),tpielemt(:,kkm,1),ichain(kkm),&      !CIP,  exclude intramolecular interaction
                            &tpicrd(:,:,kkn,2),tpielemt(:,kkn,2),ichain(kkn),&
                            & morU,morM,morV,morP,morE,morC,tpichg(:,kkm,1),tpichg(:,kkn,2),nsc,npcs,myid,ierr)
                    if (myid .eq. 0) then 
                        write(600,902,advance='no') morC
                        write(601,902,advance='no') morU
                        write(602,902,advance='no') morM
                        write(603,902,advance='no') morV
                        write(604,902,advance='no') morP
                        write(605,902,advance='no') morE
                    end if 
                end do
            end do
            if (myid .eq. 0) then 
                write(600,*);write(601,*);write(602,*);write(603,*);write(604,*);write(605,*)
            end if 
        end if ! methods pi 
    end if !method all
!end if ! frc
end do !loop

900  format(' ',i5,i5,*(xf9.5))
901  format(i5,i5,12(xf13.5),*(xf11.5))
902  format(12(xf13.5),*(xf11.5))

910 if (myid .eq. 0) write(*,*) "read connection reach end"
if (myid .eq. 0) then 
close(102)
    if (trim(method) == 'geo') then
        close(307)
        close(308)
        close(309)
        close(310)
        close(407)
        close(408)
        close(409)
        close(410)
        close(411)
    else if (trim(method) == 'geopi') then
        close(306)
        close(406)
        close(606)
        close(609)
        close(610)
    else if (trim(method) == '3dm') then
        close(300)
        close(301)
        close(302)
        close(303)
        close(304)
        close(305)
    else if (trim(method) == 'd3d') then
        close(400)
        close(401)
        close(402)
        close(403)
        close(404)
        close(405)
    else if (trim(method) == '3dpi') then
        close(500)
        close(501)
        close(502)
        close(503)
        close(504)
        close(505)
    else if (trim(method) == 'd3dpi') then
        close(600)
        close(601)
        close(602)
        close(603)
        close(604)
        close(605)
    end if 
end if 
call MPI_Finalize(ierr)
end program main
