module lib
    implicit none
contains

function multimv(a, b,tr) !!! matrix mutipile vector
    implicit none
    real, intent(in) :: a(3,3), b(3)
    real:: multimv(3)
    integer:: i,j
    character(len=1), intent(in) :: tr
    multimv=0
    if (tr == "T") then 
      do i=1,3
        do j=1,3
          multimv(i)=a(i,j)*b(j)+multimv(i)
        end do
      end do
    else 
      do i=1,3
        do j=1,3
          multimv(j)=a(i,j)*b(j)+multimv(j)
        end do
      end do
    end if 
    return
end function multimv

function cross(a,b)
    implicit none
    real, intent(in) :: a(3), b(3)
    real:: cross(3)
    cross(1)=a(2)*b(3)-a(3)*b(2) 
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
    return 
end function cross 

function dot(a,b)
    implicit none
    real, intent(in) :: a(3), b(3)
    real::dot
    dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    return 
end function dot

function angle(a,b)
    implicit none
    real, intent(in) :: a(3), b(3)
    real::angle
    angle=acos(  (a(1)*b(1)+a(2)*b(2)+a(3)*b(3)) / (sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))* &
    & sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))) )  
    return 
end function angle

function dnearestatom(a,na,b,nb)
    implicit none
    integer,intent(in)::na,nb
    real, intent(in) :: a(3,na), b(3,nb)
    real::dnearestatom,tmp(3),dis
    integer::i,j,k
    dnearestatom=100
    do i=1,na
		do j=1,nb
			tmp=a(:,i)-b(:,j)
			dis=tmp(1)*tmp(1)+tmp(2)*tmp(2)+tmp(3)*tmp(3)
			if (dis<dnearestatom) dnearestatom=dis
		end do 
    end do
    dnearestatom=sqrt(dnearestatom)
    return 
end function dnearestatom

function discen(a,na,b,nb)
    implicit none
    integer,intent(in)::na,nb
    real, intent(in) :: a(3,na), b(3,nb)
    integer::i,j,k
    real::discen,tmpa(3),tmpb(3)
    
    tmpa=0.d0;tmpb=0.d0
    do i=1,na
		tmpa=tmpa+a(:,i)
    end do
    tmpa=tmpa/na
    
    do i=1,nb
		tmpb=tmpb+b(:,i)
    end do
    tmpb=tmpb/na
    
    tmpa=tmpa-tmpb
	discen=sqrt(tmpa(1)*tmpa(1)+tmpa(2)*tmpa(2)+tmpa(3)*tmpa(3))
    return 
end function discen

end module lib 


SUBROUTINE prepareWeight(weightType, Atoms, nAtoms,weights)
implicit none 
integer:: nAtoms,i,j
character::Atoms(nAtoms)
real::weights(nAtoms)
character(2)::weightType
real::atomicMasses(29),atomicVolumes(29),atomicPolarizabilities(29),atomicElectronegativities(29)
character(2)::atomicSymbols(29)
! Data is taken from CRC Handbook of Chemistry and Physics by D.R. Lide (editor), press 2009-2010, edition. (consistent with DRAGON 6)

atomicSymbols= (/" H","Li"," C"," N"," O"," F","Na","Mg","Si"," P"," S","Cl"," K","Ni","Cu",&
& "Zn","Ga","As","Se","Br","Ag","Cd","In","Sn","Te"," I","Hg","Tl","Pb"/)
atomicMasses = (/1.01,6.941,12.01,14.01,16.0,19.0,22.991,24.305,28.09,30.97,32.07,35.45,39.098,58.69,&
& 63.55,65.39,69.72,74.92,78.96,79.9,107.87,112.41,114.82,118.71,127.6,126.9,200.59,204.38,207.2/)
atomicVolumes = (/5.42,25.25,20.58,15.6,14.71,13.31,49.0,21.69,38.79,24.43,24.43,22.45,87.11,18.14,&
& 11.49,11.25,27.39,26.52,28.73,26.52,21.31,16.52,30.11,42.8,36.62,32.52,15.6,31.54,34.53/)
atomicPolarizabilities = (/0.67,24.3,1.76,1.1,0.8,0.56,23.6,10.6,5.38,3.63,2.9,2.18,43.4,6.8,6.1,&
& 7.1,8.12,4.31,3.77,3.05,7.2,7.2,10.2,7.7,5.5,5.35,5.7,7.6,6.8/)
atomicElectronegativities =(/2.59,0.89,2.75,3.19,3.65,4.0,0.56,1.32,2.14,2.52,2.96,3.48,0.45,1.94,&
& 1.98,2.23,2.42,2.82,3.01,3.22,1.83,1.98,2.14,2.3,2.62,2.78,2.2,2.25,2.29/)
weights = 0.d0

do i =1,nAtoms
	do j = 1,29
		if (trim(adjustl(Atoms(i)))==trim(adjustl(atomicSymbols(j)))) then
			if (weightType=="ms") weights(i)=atomicMasses(j)/atomicMasses(2)!"mass":
			if (weightType=="vl") weights(i)=atomicVolumes(j)/atomicVolumes(2) !"volume":
			if (weightType=="pl") weights(i)=atomicPolarizabilities(j)/atomicPolarizabilities(2)!polarizability
			if (weightType=="en") weights(i)=atomicElectronegativities(j)/atomicElectronegativities(2) !electronegativity
			exit
		end if
	end do
end do
end SUBROUTINE


SUBROUTINE morse3d(sc,method,c,elemt,nam,morU,morM,morV,morP,morE,chg,morC)
use lib 
implicit none 
integer,intent(in)::nam
real,intent(out)::morU(32),morM(32),morV(32),morP(32),morE(32)
real,intent(in)::c(3,nam),sc
character(2),intent(in)::elemt(nam) 
character(15),intent(in)::method 
real,intent(in),optional::chg(1000)
real,intent(out),optional::morC(32)
integer::i,j,k
real::tc(3),l,r
character(10)::flname,mname,tmname
character(5)::telemnt
real::weightM(1000),weightV(1000),weightP(1000),weightE(1000),member,cwm,cwv,cwp,cwe,cwc,ttt

morU=0.d0;morM=0.d0;morV=0.d0;morP=0.d0;morE=0.d0
if (present(chg)) morC=0.d0
!3D MORSE WEIGHT calculation 
call prepareWeight("ms",elemt(1:nam),nam,weightM)
call prepareWeight("vl",elemt(1:nam),nam,weightV)
call prepareWeight("pl",elemt(1:nam),nam,weightP)
call prepareWeight("en",elemt(1:nam),nam,weightE)
L=10/3.1415926 !nm, box size

do i =2,nam
	do j=1,i-1
		cwm=weightM(i)*weightM(j)
		cwv=weightV(i)*weightV(j)
		cwp=weightP(i)*weightP(j)
		cwe=weightE(i)*weightE(j)
		if (present(chg)) cwc=chg(i)*chg(j)
		tc(1:3)=c(1:3,j)-c(1:3,i)
		r=sqrt(dot(tc,tc)) ! from nm to A
		
		do k =0,31
			if(k==0) then 
				member=1.0; 
			else
				ttt=k*r*sc
				if (trim(method)=="snsr") then  ! sin
					member=sin(ttt)/ttt
				else if (trim(method)=="snr") then
					member=sin(ttt)/r
				else if (trim(method)=="snsqrtr") then
					member=sin(ttt)/sqrt(r)
				else if (trim(method)=="sn") then
					member=sin(ttt)
				else if (trim(method)=="cnsr") then  !cos
					member=cos(ttt)/ttt
				else if (trim(method)=="cnr") then
					member=cos(ttt)/r
				else if (trim(method)=="cnsqrtr") then
					member=cos(ttt)/sqrt(r)
				else if (trim(method)=="cn") then
					member=cos(ttt)
				else if (trim(method)=="abssnsr") then  ! abssin
					member=abs(sin(ttt))/ttt
				else if (trim(method)=="abssnr") then
					member=abs(sin(ttt))/r
				else if (trim(method)=="abssnsqrtr") then
					member=abs(sin(ttt))/sqrt(r)
				else if (trim(method)=="abssn") then
					member=abs(sin(ttt))
				else if (trim(method)=="abscnsr") then  !abscos
					member=abs(cos(ttt))/ttt
				else if (trim(method)=="abscnr") then
					member=abs(cos(ttt))/r
				else if (trim(method)=="abscnsqrtr") then
					member=abs(cos(ttt))/sqrt(r)
				else if (trim(method)=="abscn") then
					member=abs(cos(ttt))
				else 
					write(*,*) "3DMoRSE core function test. Please specify the core function:"
					write(*,*) "Options:"
					write(*,*) "snsr snr sn snsqrtr and correponding cosine and abs"
					write(*,*) "As the value you given is not support, use default snsr"
					member=sin(ttt)/ttt
				end if 
			end if 
			morU(k+1)=morU(k+1)+member
			morM(k+1)=morM(k+1)+member*cwm
			morV(k+1)=morV(k+1)+member*cwv
			morP(k+1)=morP(k+1)+member*cwp
			morE(k+1)=morE(k+1)+member*cwe
			if (present(chg)) morC(k+1)=morC(k+1)+member*cwc
		end do
	end do
end do 
END SUBROUTINE 


!!!! only for interaction remove intramolecular interaction 
!!!!CIP !!! 
SUBROUTINE morse3ddimer(sc,method,c1,elemt1,nam1,c2,elemt2,nam2,morU,morM,morV,morP,morE,chg1,chg2,morC,cip)
use lib 
implicit none 
integer,intent(in)::nam1,nam2
real,intent(in)::c1(3,nam1),c2(3,nam2),sc
character(15),intent(in)::method 
character(2),intent(in)::elemt1(nam1),elemt2(nam2)
real,intent(in),optional::chg1(1000),chg2(1000)
real,intent(out),optional::morC(32)
real,intent(out)::morU(32),morM(32),morV(32),morP(32),morE(32),cip(nam1,nam2)
integer::i,j,k
real::tc(3),l,r
character(10)::flname,mname,tmname
character(5)::telemnt
real::weightM1(1000),weightV1(1000),weightP1(1000),weightE1(1000),member,cwm,cwv,cwp,cwe,cwc,ttt
real::weightM2(1000),weightV2(1000),weightP2(1000),weightE2(1000)

morU=0.d0;morM=0.d0;morV=0.d0;morP=0.d0;morE=0.d0;cip=0.d0
if (present(chg1)) morC=0.d0

!3D MORSE WEIGHT calculation 
call prepareWeight("ms",elemt1(1:nam1),nam1,weightM1)
call prepareWeight("vl",elemt1(1:nam1),nam1,weightV1)
call prepareWeight("pl",elemt1(1:nam1),nam1,weightP1)
call prepareWeight("en",elemt1(1:nam1),nam1,weightE1)

call prepareWeight("ms",elemt2(1:nam2),nam2,weightM2)
call prepareWeight("vl",elemt2(1:nam2),nam2,weightV2)
call prepareWeight("pl",elemt2(1:nam2),nam2,weightP2)
call prepareWeight("en",elemt2(1:nam2),nam2,weightE2)

do i =1,nam1
	do j=1,nam2
		cwm=weightM1(i)*weightM2(j)
		cwv=weightV1(i)*weightV2(j)
		cwp=weightP1(i)*weightP2(j)
		cwe=weightE1(i)*weightE2(j)
		if (present(chg1)) cwc=chg1(i)*chg2(j)
		
		tc(1:3)=c2(1:3,j)-c1(1:3,i)
		r=sqrt(dot(tc,tc)) ! from nm to A
		
		if (present(chg1)) cip(i,j)=cwc/r  !! cip CIP
		
		do k =0,31
			if(k==0) then 
				member=1.0; 
			else
				ttt=k*r*sc
				if (trim(method)=="snsr") then  ! sin
					member=sin(ttt)/ttt
				else if (trim(method)=="snr") then
					member=sin(ttt)/r
				else if (trim(method)=="snsqrtr") then
					member=sin(ttt)/sqrt(r)
				else if (trim(method)=="sn") then
					member=sin(ttt)
				else if (trim(method)=="cnsr") then  !cos
					member=cos(ttt)/ttt
				else if (trim(method)=="cnr") then
					member=cos(ttt)/r
				else if (trim(method)=="cnsqrtr") then
					member=cos(ttt)/sqrt(r)
				else if (trim(method)=="cn") then
					member=cos(ttt)
				else if (trim(method)=="abssnsr") then  ! abssin
					member=abs(sin(ttt))/ttt
				else if (trim(method)=="abssnr") then
					member=abs(sin(ttt))/r
				else if (trim(method)=="abssnsqrtr") then
					member=abs(sin(ttt))/sqrt(r)
				else if (trim(method)=="abssn") then
					member=abs(sin(ttt))
				else if (trim(method)=="abscnsr") then  !abscos
					member=abs(cos(ttt))/ttt
				else if (trim(method)=="abscnr") then
					member=abs(cos(ttt))/r
				else if (trim(method)=="abscnsqrtr") then
					member=abs(cos(ttt))/sqrt(r)
				else if (trim(method)=="abscn") then
					member=abs(cos(ttt))
				else 
					write(*,*) "3DMoRSE core function test. Please specify the core function:"
					write(*,*) "Options:"
					write(*,*) "snsr snr sn snsqrtr and correponding cosine and abs"
					write(*,*) "As the value you given is not support, use default snsr"
					member=sin(ttt)/ttt
				end if 
			end if 
			morU(k+1)=morU(k+1)+member
			morM(k+1)=morM(k+1)+member*cwm
			morV(k+1)=morV(k+1)+member*cwv
			morP(k+1)=morP(k+1)+member*cwp
			morE(k+1)=morE(k+1)+member*cwe
			if (present(chg1)) morC(k+1)=morC(k+1)+member*cwc
		end do
	end do
end do 
!write(*,900) cip 
!read(*,*) i 
!900  format(' ',*(xf13.5))
END SUBROUTINE 

! 3D-MORSE based on nearest distance of each atoms
SUBROUTINE msnearest(c1,elemt1,nam1,c2,elemt2,nam2,morU,morM,morV,morP,morE,chg1,chg2,morC)
use lib 
implicit none 
integer,intent(in)::nam1,nam2
real,intent(in)::c1(3,nam1),c2(3,nam2)
character(2),intent(in)::elemt1(nam1),elemt2(nam2) 
real,intent(in),optional::chg1(1000),chg2(1000)
real,intent(out),optional::morC(32)
real,intent(out)::morU(32),morM(32),morV(32),morP(32),morE(32)
integer::i,j,k
real::tc(3),l,r,rtmp,itmp
character(10)::flname,mname,tmname
character(5)::telemnt
real::weightM1(1000),weightV1(1000),weightP1(1000),weightE1(1000),member,cwm,cwv,cwp,cwe,cwc,ttt
real::weightM2(1000),weightV2(1000),weightP2(1000),weightE2(1000)

morU=0.d0;morM=0.d0;morV=0.d0;morP=0.d0;morE=0.d0
if (present(chg1)) morC=0.d0

!3D MORSE WEIGHT calculation 
call prepareWeight("ms",elemt1(1:nam1),nam1,weightM1)
call prepareWeight("vl",elemt1(1:nam1),nam1,weightV1)
call prepareWeight("pl",elemt1(1:nam1),nam1,weightP1)
call prepareWeight("en",elemt1(1:nam1),nam1,weightE1)

call prepareWeight("ms",elemt2(1:nam2),nam2,weightM2)
call prepareWeight("vl",elemt2(1:nam2),nam2,weightV2)
call prepareWeight("pl",elemt2(1:nam2),nam2,weightP2)
call prepareWeight("en",elemt2(1:nam2),nam2,weightE2)

do i =1,nam1
	rtmp=100.0
	itmp=1
	do j=1,nam2
		tc(1:3)=c2(1:3,j)-c1(1:3,i)
		r=sqrt(dot(tc,tc)) ! from nm to A
		if (rtmp> r) then 
			rtmp=r ! find nearest
			itmp=j
		end if 
	end do
	j=itmp
	r=rtmp 
	cwm=weightM1(i)*weightM2(j)
	cwv=weightV1(i)*weightV2(j)
	cwp=weightP1(i)*weightP2(j)
	cwe=weightE1(i)*weightE2(j)
	if (present(chg1)) cwc=chg1(i)*chg2(j)
	do k =0,31
		if(k==0) then 
			member=1.0
		else
			ttt=k*r*10.0
			member=sin(ttt)/ttt
			!member=member*member
		end if 
		morU(k+1)=morU(k+1)+member
		morM(k+1)=morM(k+1)+member*cwm
		morV(k+1)=morV(k+1)+member*cwv
		morP(k+1)=morP(k+1)+member*cwp
		morE(k+1)=morE(k+1)+member*cwe
		if (present(chg1)) morC(k+1)=morC(k+1)+member*cwc
	end do
end do 
END SUBROUTINE 

! 32th of nearest distance 
SUBROUTINE dnearest(c1,elemt1,nam1,c2,elemt2,nam2,morU,morM,morV,morP,morE,chg1,chg2,morC)
use lib 
implicit none 
integer,intent(in)::nam1,nam2
real,intent(in)::c1(3,nam1),c2(3,nam2)
character(2),intent(in)::elemt1(nam1),elemt2(nam2) 
real,intent(in),optional::chg1(1000),chg2(1000)
real,intent(out),optional::morC(32)
real,intent(out)::morU(32),morM(32),morV(32),morP(32),morE(32)
integer::i,j,k
real::tc(3),l,r,rtmp,itmp
character(10)::flname,mname,tmname
character(5)::telemnt
real::weightM1(1000),weightV1(1000),weightP1(1000),weightE1(1000),member,cwm,cwv,cwp,cwe,cwc,ttt
real::weightM2(1000),weightV2(1000),weightP2(1000),weightE2(1000)

morU=100.d0;morM=0.d0;morV=0.d0;morP=0.d0;morE=0.d0
if (present(chg1)) morC=0.d0

!3D MORSE WEIGHT calculation 
call prepareWeight("ms",elemt1(1:nam1),nam1,weightM1)
call prepareWeight("vl",elemt1(1:nam1),nam1,weightV1)
call prepareWeight("pl",elemt1(1:nam1),nam1,weightP1)
call prepareWeight("en",elemt1(1:nam1),nam1,weightE1)

call prepareWeight("ms",elemt2(1:nam2),nam2,weightM2)
call prepareWeight("vl",elemt2(1:nam2),nam2,weightV2)
call prepareWeight("pl",elemt2(1:nam2),nam2,weightP2)
call prepareWeight("en",elemt2(1:nam2),nam2,weightE2)

do i =1,nam1
	do j=1,nam2
		tc(1:3)=c2(1:3,j)-c1(1:3,i)
		r=sqrt(dot(tc,tc)) ! from nm to A
		do k =1,32
			if (morU(k)> r) then 
				if (k>1) morU(k-1)=morU(k)
				morU(k)=r
			end if 
		end do 
	end do
end do 
END SUBROUTINE 

! distance of 4v, output 16 distance 
SUBROUTINE  calparam4(indx1,crd1,natom1,indx2,crd2,natom2,dis)
use lib 
implicit none
integer,intent(in)::indx1(4),natom1,indx2(4),natom2
integer::i,j,k
real::crd1(3,natom1),crd2(3,natom2),vtmp(3)
real,intent(out)::dis(16)
k=1
do i =1,4
	do j=1,4
	    vtmp(:)=crd1(:,indx1(i))-crd2(:,indx2(j))
	    dis(k)=dot(vtmp,vtmp)
	    k=k+1
	end do
end do 
END SUBROUTINE 

!parameters based on 4 chosen points
SUBROUTINE  calparam4v(indx1,crd1,natom1,indx2,crd2,natom2,dis)
use lib 
implicit none
integer,intent(in)::indx1(4),natom1,indx2(4),natom2
integer::i,j,k
real::crd1(3,natom1),crd2(3,natom2),v1(3,4),v2(3,4),vn1(3,2),vn2(3,2)
real::cen1(3),cen2(3),outparam1(10),outparam2(10),outparam(10)
real,intent(out)::dis(32)

!!! molecule center
cen1(:)=0.d0
do i=1,natom1
	cen1(:)=cen1(:)+crd1(:,i)
end do
cen1(:)=cen1(:)/natom1
cen2(:)=0.d0
do i=1,natom2
	cen2(:)=cen2(:)+crd2(:,i)
end do
cen2(:)=cen2(:)/natom2

vn1(:,1)=cen2-cen1  ! center-center vector
outparam(1)=sqrt(dot(vn1(:,1),vn1(:,1))) ! dcen
outparam(2)=dnearestatom(crd1(1:3,1:natom1),natom1,crd2(1:3,1:natom2),natom2) !datom

call calparam1r(indx1,crd1,natom1,outparam1)
call calparam1r(indx2,crd2,natom2,outparam2)

v1(:,1)=crd1(:,indx1(3))-crd1(:,indx1(1))
v1(:,2)=crd1(:,indx1(4))-crd1(:,indx1(1))
v1(:,3)=crd1(:,indx1(3))-crd1(:,indx1(2))
v1(:,4)=crd1(:,indx1(4))-crd1(:,indx1(2))

v2(:,1)=crd2(:,indx2(3))-crd2(:,indx2(1))
v2(:,2)=crd2(:,indx2(4))-crd2(:,indx2(1))
v2(:,3)=crd2(:,indx2(3))-crd2(:,indx2(2))
v2(:,4)=crd2(:,indx2(4))-crd2(:,indx2(2))

outparam(3)=dot(vn1(:,1),cross(v1(:,1), v1(:,2)))/outparam1(5)   ! dcen projection on 1,param1(:,5) normal V
outparam(4)=dot(vn1(:,1),cross(v2(:,1), v2(:,2)))/outparam2(5)   ! dcen projection on 2,param2(:,5) normal V
outparam(5)=acos(dot(v1(:,1),v2(:,1))/(outparam1(1)*outparam2(1))) !cos angle of long axis
outparam(6)=acos(dot(v1(:,2),v2(:,2))/(outparam1(2)*outparam2(2))) !cos angle of short axis

vn1(:,1)=cross(v1(:,3),v1(:,1)) !# normal of right plane 321
vn1(:,1)=vn1(:,1)/sqrt(dot(vn1(:,1),vn1(:,1)))
vn1(:,2)=cross(v1(:,4),v1(:,2)) !# normal of left plane  421
vn1(:,2)=vn1(:,2)/sqrt(dot(vn1(:,2),vn1(:,2)))

vn2(:,1)=cross(v2(:,3),v2(:,1)) !# normal of right plane 321
vn2(:,1)=vn2(:,1)/sqrt(dot(vn2(:,1),vn2(:,1)))
vn2(:,2)=cross(v2(:,4),v2(:,2)) !# normal of left plane  421
vn2(:,2)=vn2(:,2)/sqrt(dot(vn2(:,2),vn2(:,2)))

outparam(7) =acos(dot(vn1(:,1),vn2(:,1))) !cos angle of 11-21 normal axis
outparam(8) =acos(dot(vn1(:,2),vn2(:,2))) !cos angle of 12-22 normal axis
outparam(9) =acos(dot(vn1(:,1),vn2(:,2))) !cos angle of 11-22 normal axis
outparam(10)=acos(dot(vn1(:,2),vn2(:,1))) !cos angle of 12-21normal axis

dis(1:6)=outparam1(1:6)
dis(7:12)=outparam2(1:6)
dis(13:22)=outparam(1:10)
END SUBROUTINE 


SUBROUTINE  calparam2(indx1,crd1,natom1,indx2,crd2,natom2,outparamall)
use lib 
implicit none
integer,intent(in)::indx1(4),natom1,indx2(4),natom2
integer::i
real::crd1(3,natom1),crd2(3,natom2),param1(3,5),vvtmp(3,2),param2(3,5)
real::outparam1(10),outparam2(10),outparam(10)
real,intent(out)::outparamall(32)

!!!!!!!!!!!!!!!!!!!!molecule 1 !!!!!!!!!!!!!!!!!!!!!!!!!
param1(:,1)=crd1(:,indx1(2))-crd1(:,indx1(1)) !# long axis
param1(:,2)=crd1(:,indx1(4))-crd1(:,indx1(3)) ! short
param1(:,3)=crd1(:,indx1(3))-crd1(:,indx1(1)) !# right
param1(:,4)=crd1(:,indx1(4))-crd1(:,indx1(1)) !# left
param1(:,5)=cross(param1(:,1), param1(:,2)) !long X short = normal,area
call calparam1(indx1,crd1,natom1,outparam1)
!!!!!!!!!!!!!!!!!!!!molecule 2 !!!!!!!!!!!!!!!!!!!!!!!!!
param2(:,1)=crd2(:,indx2(2))-crd2(:,indx2(1)) !# long axis
param2(:,2)=crd2(:,indx2(4))-crd2(:,indx2(3)) ! short
param2(:,3)=crd2(:,indx2(3))-crd2(:,indx2(1)) !# right
param2(:,4)=crd2(:,indx2(4))-crd2(:,indx2(1)) !# left
param2(:,5)=cross(param2(:,1), param2(:,2)) !long X short = normal,area
call calparam1(indx2,crd2,natom2,outparam2)
!write(*,*) "output3",outparam2(3:5)
!!!!!!!!!!dimer !!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!! molecule center
vvtmp(:,1)=0.d0
do i=1,natom1
	vvtmp(:,1)=vvtmp(:,1)+crd1(:,i)
end do
vvtmp(:,1)=vvtmp(:,1)/natom1

vvtmp(:,2)=0.d0
do i=1,natom2
	vvtmp(:,2)=vvtmp(:,2)+crd2(:,i)
end do
vvtmp(:,2)=vvtmp(:,2)/natom2

vvtmp(:,1)=vvtmp(:,2)-vvtmp(:,1)  ! center-center vector
outparam(1)=sqrt(dot(vvtmp(:,1),vvtmp(:,1))) ! dcen
outparam(2)=dnearestatom(crd1(1:3,1:natom1),natom1,crd2(1:3,1:natom2),natom2) !datom
outparam(3)=dot(vvtmp(:,1),param1(:,5))/sqrt(dot(param1(:,5),param1(:,5)) )  ! dcen projection on 1,param1(:,5) normal V
outparam(4)=dot(vvtmp(:,1),param2(:,5))/sqrt(dot(param2(:,5),param2(:,5)) )   ! dcen projection on 2,param2(:,5) normal V

outparam(5)=dot(param1(:,1),param2(:,1))/sqrt(dot(param1(:,1),param1(:,1))*dot(param2(:,1),param2(:,1)))   !cos angle of long axis
outparam(6)=dot(param1(:,2),param2(:,2))/sqrt(dot(param1(:,2),param1(:,2))*dot(param2(:,2),param2(:,2)))  !cos angle of short axis
outparam(7)=dot(param1(:,3),param2(:,3))/sqrt(dot(param1(:,5),param1(:,5))*dot(param2(:,5),param2(:,5)))   !cos angle of normal axis

vvtmp(:,1)=cross(param1(:,1),param1(:,3))   ! orthogonal of long and normal axis
vvtmp(:,2)=cross(param2(:,1),param2(:,3))  
vvtmp(:,1)=vvtmp(:,1)/sqrt(dot(vvtmp(:,1),vvtmp(:,1))) ! nomalize
vvtmp(:,2)=vvtmp(:,2)/sqrt(dot(vvtmp(:,2),vvtmp(:,2)))   
outparam(8)=dot(vvtmp(:,1),vvtmp(:,2))        !cos angle 

vvtmp(:,1)= dot(param1(:,1),param2(:,5))*param2(:,5)/dot(param2(:,5),param2(:,5)) ! L1 on N2, P12 always has same direction to L1
vvtmp(:,2)= param1(:,1) - vvtmp(:,1)  ! L1 projection on plane 
!            S1          *S2         * COS normal  * COS projection of L
outparam(9)=sqrt(dot(param1(:,5),param1(:,5)))*sqrt(dot(param2(:,5),param2(:,5)))* &
& outparam(7)  *dot(param2(:,1),vvtmp(:,2))/sqrt(dot(param2(:,2),param2(:,2)))/sqrt(dot(vvtmp(:,2),vvtmp(:,2)))

outparamall(1:6)=outparam1(1:6)
outparamall(7:12)=outparam2(1:6)
outparamall(13:21)=outparam(1:9)
END SUBROUTINE 


SUBROUTINE  calparam2r(indx1,crd1,natom1,indx2,crd2,natom2,outparamall)  !use param1 for molecule,not good
use lib 
implicit none
integer,intent(in)::indx1(4),natom1,indx2(4),natom2
integer::i
real::crd1(3,natom1),crd2(3,natom2),param1(3,5),vvtmp(3,2),param2(3,5)
real::outparam1(10),outparam2(10),outparam(10)
real,intent(out)::outparamall(32)
!!!!!!!!!!!!!!!!!!!!molecule 1 !!!!!!!!!!!!!!!!!!!!!!!!!
param1(:,1)=crd1(:,indx1(2))-crd1(:,indx1(1)) !# long axis
param1(:,2)=crd1(:,indx1(4))-crd1(:,indx1(3)) ! short
param1(:,3)=crd1(:,indx1(3))-crd1(:,indx1(1)) !# right
param1(:,4)=crd1(:,indx1(4))-crd1(:,indx1(1)) !# left
param1(:,5)=cross(param1(:,1), param1(:,2)) !long X short = normal,area
call calparam1r(indx1,crd1,natom1,outparam1)
!!!!!!!!!!!!!!!!!!!!molecule 2 !!!!!!!!!!!!!!!!!!!!!!!!!
param2(:,1)=crd2(:,indx2(2))-crd2(:,indx2(1)) !# long axis
param2(:,2)=crd2(:,indx2(4))-crd2(:,indx2(3)) ! short
param2(:,3)=crd2(:,indx2(3))-crd2(:,indx2(1)) !# right
param2(:,4)=crd2(:,indx2(4))-crd2(:,indx2(1)) !# left
param2(:,5)=cross(param2(:,1), param2(:,2)) !long X short = normal,area
call calparam1r(indx2,crd2,natom2,outparam2)
!write(*,*) "output3",outparam2(3:5)
!!!!!!!!!!dimer !!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!! molecule center
vvtmp(:,1)=0.d0
do i=1,natom1
	vvtmp(:,1)=vvtmp(:,1)+crd1(:,i)
end do
vvtmp(:,1)=vvtmp(:,1)/natom1

vvtmp(:,2)=0.d0
do i=1,natom2
	vvtmp(:,2)=vvtmp(:,2)+crd2(:,i)
end do
vvtmp(:,2)=vvtmp(:,2)/natom2

vvtmp(:,1)=vvtmp(:,2)-vvtmp(:,1)  ! center-center vector
outparam(1)=sqrt(dot(vvtmp(:,1),vvtmp(:,1))) ! dcen
outparam(2)=dnearestatom(crd1(1:3,1:natom1),natom1,crd2(1:3,1:natom2),natom2) !datom
outparam(3)=dot(vvtmp(:,1),param1(:,5))/sqrt(dot(param1(:,5),param1(:,5)) )  ! dcen projection on 1,param1(:,5) normal V
outparam(4)=dot(vvtmp(:,1),param2(:,5))/sqrt(dot(param2(:,5),param2(:,5)) )   ! dcen projection on 2,param2(:,5) normal V

outparam(5)=dot(param1(:,1),param2(:,1))/sqrt(dot(param1(:,1),param1(:,1))*dot(param2(:,1),param2(:,1)))   !cos angle of long axis
outparam(6)=dot(param1(:,2),param2(:,2))/sqrt(dot(param1(:,2),param1(:,2))*dot(param2(:,2),param2(:,2)))  !cos angle of short axis
outparam(7)=dot(param1(:,3),param2(:,3))/sqrt(dot(param1(:,5),param1(:,5))*dot(param2(:,5),param2(:,5)))   !cos angle of normal axis

vvtmp(:,1)=cross(param1(:,1),param1(:,3))   ! orthogonal of long and normal axis
vvtmp(:,2)=cross(param2(:,1),param2(:,3))  
vvtmp(:,1)=vvtmp(:,1)/sqrt(dot(vvtmp(:,1),vvtmp(:,1))) ! nomalize
vvtmp(:,2)=vvtmp(:,2)/sqrt(dot(vvtmp(:,2),vvtmp(:,2)))   
outparam(8)=dot(vvtmp(:,1),vvtmp(:,2))        !cos angle 

vvtmp(:,1)= dot(param1(:,1),param2(:,5))*param2(:,5)/dot(param2(:,5),param2(:,5)) ! L1 on N2, P12 always has same direction to L1
vvtmp(:,2)= param1(:,1) - vvtmp(:,1)  ! L1 projection on plane 
!            S1          *S2         * COS normal  * COS projection of L
outparam(9)=sqrt(dot(param1(:,5),param1(:,5)))*sqrt(dot(param2(:,5),param2(:,5)))* &
& outparam(7)  *dot(param2(:,1),vvtmp(:,2))/sqrt(dot(param2(:,2),param2(:,2)))/sqrt(dot(vvtmp(:,2),vvtmp(:,2)))

outparamall(1:6)=outparam1(1:6)
outparamall(7:12)=outparam2(1:6)
outparamall(13:21)=outparam(1:9)
END SUBROUTINE 

!! Lederer's discriptors
SUBROUTINE  calref(indx1,natom1,crd1,elemt1,indx2,natom2,crd2,elemt2,outparam)
use lib 
implicit none
integer,intent(in)::indx1(4),natom1,indx2(4),natom2
character(2),intent(in)::elemt1(natom1),elemt2(natom2)
integer::i
real::crd1(3,natom1),crd2(3,natom2),mcen1(3),mcen2(3),tm,vtmp(3)
real::weightM(1000),param2(3,5),param1(3,5),outparam(6)

!!!!!!!!!!!!!!!!!!!!molecule 1 !!!!!!!!!!!!!!!!!!!!!!!!!
call prepareWeight("ms",elemt1(1:natom1),natom1,weightM)
mcen1=0
tm=0
do i=1,natom1
	mcen1(1:3)=mcen1(1:3)+crd1(1:3,i)*weightM(i)
	tm=tm+weightM(i)
end do
mcen1(1:3)=mcen1(1:3)/tm

call prepareWeight("ms",elemt2(1:natom2),natom2,weightM)
mcen2=0
tm=0
do i=1,natom2
	mcen2(1:3)=mcen2(1:3)+crd2(1:3,i)*weightM(i)
	tm=tm+weightM(i)
end do
mcen2(1:3)=mcen2(1:3)/tm

mcen1=mcen2-mcen1  ! r
param1(:,1)=crd1(:,indx1(2))-crd1(:,indx1(1)) !# long axis
param1(:,2)=crd1(:,indx1(4))-crd1(:,indx1(3)) ! short
param1(:,3)=cross(param1(:,1), param1(:,2)) !long X short = normal,area
param1(:,2)=cross(param1(:,3),param1(:,1))   ! short, same with reference
!nomalized
param1(:,1)=param1(:,1)/sqrt(dot(param1(:,1),param1(:,1)))
param1(:,2)=param1(:,2)/sqrt(dot(param1(:,2),param1(:,2)))
param1(:,3)=param1(:,3)/sqrt(dot(param1(:,3),param1(:,3)))
!position
vtmp(1)=dot(mcen1,param1(:,2)) ! x- short prejection
vtmp(2)=dot(mcen1,param1(:,3)) ! y- normal
if (vtmp(2)<0) then 
	vtmp(2)=-vtmp(2)  ! it must be >0 
	vtmp(1)=-vtmp(1)  ! if y change sign, keep z unchanged, x has to change too
end if
vtmp(3)=dot(mcen1,param1(:,1)) ! z- long
if (vtmp(3)<0) then 
	vtmp(3)=-vtmp(3)  ! it must be >0 
	vtmp(1)=-vtmp(1)  ! if z change sign, keep y unchanged, x has to change too
end if
outparam(1:3)=vtmp(1:3)
!short
!!!!!!!!!!!!!!!!!!!!molecule 2 !!!!!!!!!!!!!!!!!!!!!!!!!
param2(:,1)=crd2(:,indx2(2))-crd2(:,indx2(1)) !# long axis
param2(:,2)=crd2(:,indx2(4))-crd2(:,indx2(3)) ! short
param2(:,3)=cross(param2(:,1),param2(:,2))    ! long X short = normal,area
param2(:,2)=cross(param2(:,3),param2(:,1))    ! short, same with reference
!nomalized
param2(:,1)=param2(:,1)/sqrt(dot(param2(:,1),param2(:,1)))
param2(:,2)=param2(:,2)/sqrt(dot(param2(:,2),param2(:,2)))
param2(:,3)=param2(:,3)/sqrt(dot(param2(:,3),param2(:,3)))

!write(*,*) param1(:,1),param1(:,2),param1(:,3)
!write(*,*) param2(:,1),param2(:,2),param2(:,3)
outparam(4)=acos(min(dot(param2(:,1),param1(:,1)),1.d0)) ! long angle 
outparam(5)=acos(min(dot(param2(:,2),param1(:,2)),1.d0)) ! short angle 
outparam(6)=acos(min(dot(param2(:,3),param1(:,3)),1.d0)) ! normal angle 
END SUBROUTINE 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! for mono molecules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calparam1(indx,crd,natom,outparam)
use lib 
implicit none
integer,intent(in)::indx(4),natom
integer:: i,j,k
real,intent(in)::crd(3,natom)
real::sumh,maxh,height,param(3,5),vvtmp(3,2),cen(3)
real,intent(out)::outparam(10)
param(:,1)=crd(:,indx(2))-crd(:,indx(1)) !# long axis
param(:,2)=crd(:,indx(4))-crd(:,indx(3)) ! short
param(:,3)=crd(:,indx(3))-crd(:,indx(1)) !# right
param(:,4)=crd(:,indx(4))-crd(:,indx(1)) !# left
param(:,5)=cross(param(:,1), param(:,2)) !long X short = normal,area
outparam(1)=sqrt(dot(param(:,1),param(:,1)))  !# length of long 
outparam(2)=sqrt(dot(param(:,2),param(:,2)))  !# length of short 
outparam(5)=sqrt(dot(param(:,5),param(:,5)))  ! area 
vvtmp(:,1)=cross(param(:,3),param(:,1)) !# normal of right plane
vvtmp(:,2)=cross(param(:,4),param(:,1)) !# normal of left plane
outparam(6)=-dot(vvtmp(:,1),vvtmp(:,2))/sqrt(dot(vvtmp(:,1),vvtmp(:,1))*dot(vvtmp(:,2),vvtmp(:,2))) !angle between left and right normal V
outparam(3)=dot(param(:,3),param(:,1))/outparam(1)  !# hight of right  A*B=|A||B|COS, H=ACOS=A*B/|B|
outparam(4)=dot(param(:,4),param(:,1))/outparam(1)  !# hight of left
outparam(3)=dot(param(:,3),param(:,5))/outparam(5)  !# hight of right  A*B=|A||B|COS, H=ACOS=A*B/|B|
outparam(4)=dot(param(:,4),param(:,5))/outparam(5)  !# hight of left
END SUBROUTINE 

SUBROUTINE calparam1r(indx,crd,natom,outparam) !not good
use lib 
implicit none
integer,intent(in)::indx(4),natom
integer:: i,j,k
real,intent(in)::crd(3,natom)
real::sumh,maxh,height,param(3,5),vvtmp(3,2),cen(3)
real,intent(out)::outparam(10)
param(:,1)=crd(:,indx(3))-crd(:,indx(1)) !# long axis
param(:,2)=crd(:,indx(4))-crd(:,indx(1)) ! short
param(:,3)=crd(:,indx(3))-crd(:,indx(2)) !# long 2
param(:,4)=crd(:,indx(4))-crd(:,indx(2)) !# short 2
param(:,5)=cross(param(:,1), param(:,2)) !long X short = normal,area
outparam(1)=sqrt(dot(param(:,1),param(:,1)))  !# length of long 
outparam(2)=sqrt(dot(param(:,2),param(:,2)))  !# length of short 
outparam(5)=sqrt(dot(param(:,5),param(:,5)))  ! area 
vvtmp(:,1)=cross(param(:,3),param(:,1)) !# normal of right plane 321
vvtmp(:,1)=vvtmp(:,1)/sqrt(dot(vvtmp(:,1),vvtmp(:,1)))
vvtmp(:,2)=cross(param(:,4),param(:,2)) !# normal of left plane  421
vvtmp(:,2)=vvtmp(:,2)/sqrt(dot(vvtmp(:,2),vvtmp(:,2)))
outparam(6)=-dot(vvtmp(:,1),vvtmp(:,2))!angle between left and right normal V
outparam(3)=dot(param(:,3),vvtmp(:,2))  !# hight of right  3 in 421 plane
outparam(4)=dot(param(:,4),vvtmp(:,1)) !# hight of left 
END SUBROUTINE 

