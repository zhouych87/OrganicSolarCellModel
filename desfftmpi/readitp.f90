
SUBROUTINE readitp(itpfile,nchain,chain,ichain,chg,nnc,nnca,natom)
!SUBROUTINE readitp(itpfile,nchain,lbchain,ichain,chg,nnc,lnca)
implicit none 
character(20),intent(in)::itpfile
integer,intent(out)::nchain,chain(1000,10),ichain(10)
real,intent(out),optional:: chg(1000)  !If (Present (chg)) Then 

integer::i,j,k,resid(1000),cgnr(1000),natom,nbond,napi,nbpi
integer::ii,iii,nnc,nnca(100),tnc
real::mass(1000)
character(2)::elemnt(1000)
character(6)::label(1000),lbchain(1000,10),lnca(100)
character(4)::res(1000),etype(1000)
character(20)::fformat
character(54)::linetxt
real,allocatable::new_chg(:)
integer,allocatable::anb(:),ab(:,:),cpi(:),pi_index(:),color(:)

open(101,file=itpfile,status='old') 

do while (linetxt(1:9) .ne. '[ atoms ]')
	read(101,"(a)") linetxt  !  
end do 
i=0

!write(*,*) " [ atoms ] .vs.",linetxt
read(101,"(a)") linetxt  !nr  type  resnr  resid  atom  cgnr  charge    mass
do 
	i=i+1            !! 1    HC    1    SIW4     H9    1    0.128   1.0080
	read(101,*,err=210) k, etype(i),resid(i),res(i),label(i),cgnr(i),chg(i),mass(i)
end do 

210 backspace(101)
!write(*,*) k, etype(i),resid(i),res(i),elemnt(i),cgnr(i),chg(i),mass(i)
read(101,"(a)") linetxt
if (linetxt(1:1) .ne. ';') then
	write(*,*) "REad atom error",i,k, etype(i),res(i),resid(i),label(i),cgnr(i),chg(i),mass(i)
	stop
else 
	write(*,*) "Read atom finished"
	write(*,*) "Going to read bond"
end if

natom=i-1
write(*,*) "Number of atoms: ", natom

allocate(anb(natom),ab(8,natom),cpi(natom),pi_index(natom),new_chg(natom))
new_chg=chg(1:natom)

read(101,"(a)") linetxt  !   [ bonds ] 
!write(*,*) " [ bonds ] .vs.",linetxt
read(101,"(a)") linetxt  !  ;  ai   aj  funct   c0         c1
!write(*,*) ";  ai   aj  funct   c0         c1 .vs. ",linetxt

anb=0
ab=0
i=0
do
	i=i+1 
	read(101,*,err=220) j,k 
	anb(j)=anb(j)+1 ! bond number/index of atom bond(1)
	anb(k)=anb(k)+1 ! bond number/index of atom bond(2)
	ab(anb(j),j)=k ! the anb bond of atom bond(1) is connected with atom-bond(2)
	ab(anb(k),k)=j ! the anb bond of atom bond(2) is connected with atom-bond(1)
end do

220 backspace(101)
read(101,"(a)") linetxt
if (linetxt(1:1) .ne. '[') then  ![ pairs ]
	write(*,*) "Read band error",i,j,k
	stop
else 
	write(*,*) "Read bond finished"
end if

nbond=i-1
close(101)

napi=0
nnc=0
do i =1,natom! ! not C, not H 
	if  ((etype(i)(1:1).ne. "C") .and. (etype(i)(1:1) .ne. "H") ) then  ! caution: Cl is not include
		cpi(i)=2             ! is part of conjugation 
		napi=napi+1
		pi_index(napi)=i    !! i in molecule -> napi in pi 
	else if ((etype(i)(1:1) .eq. "C") .and. & 
	& ((etype(i)(2:2) .ne. "R") .and. (etype(i)(2:2) .ne. " ") ) )  then !!! not CH3 CR3, thi include Cl,CH3 is C?? 
		cpi(i)=1             ! is part of conjugation 
		napi=napi+1
		pi_index(napi)=i    !! i in molecule -> napi in pi 
	else
		cpi(i)=0
	end if 
end do 

write(*,*) "number of pi atoms in the whole molecule",napi
write(*,"(a)",ADVANCE="NO") "Whether atoms are pi atoms"
do i=1,natom
	write(*,"(xi0)",ADVANCE="NO") cpi(i)
end do 
write(*,*) 
write(*,"(a)",ADVANCE="NO") "pi atoms"
do i=1,napi
	write(*,"(xi0)",ADVANCE="NO") pi_index(i)
end do 
write(*,*) 


allocate(color(natom))
color=-1
nchain=0

do iii =1,natom ! assume there are napi conjugation planes
	do i =1,natom
		if (cpi(i)>=1 .and. color(i)==-1) then 
			nchain=nchain+1
			color(i)=nchain
			exit
		end if
	end do
!!! before set initial color, below spread color to others
	do k =1,nbond !  the maximum loop should be less than the number of bonds
		do i =1,natom
			if (cpi(i)>=1 .and. color(i) .ne. -1) then
				do j =1,anb(i)
					if (cpi(ab(j,i))>=1) color(ab(j,i))=color(i)
				end do
			end if
		end do
	end do 
end do 

!end of find planes
write(*,*) "nchain: ",nchain
write(*,'(a)',ADVANCE="NO") "color: "
do i=1,natom
	write(*,"(xi0)",ADVANCE="NO") color(i)
end do 
write(*,*) 

!allocate(chain(napi-nchain,nchain),ichain(nchain))
chain=0
ichain=0

do i =1,napi
	ichain(color(pi_index(i)))=ichain(color(pi_index(i)))+1
	chain(ichain(color(pi_index(i))),color(pi_index(i)))=pi_index(i)
	lbchain(ichain(color(pi_index(i))),color(pi_index(i)))=label(pi_index(i))
end do

tnc=0
nnc=0 ! number of non-carbon atoms
do i=1,nchain
	write(*,*) "ichain: ",i,ichain(i) 
	write(*,*) "chain: ",i
	do j=1,ichain(i)
		write(*,"(xi0)",ADVANCE="NO") chain(j,i)
	end do 
	write(*,*) 
	
	if (ichain(i)>=4) then
		tnc=tnc+1
		ichain(tnc)=ichain(i) ! tnc <= i 
		chain(:,tnc)=chain(:,i)
	else 
		do j=1,ichain(i)
			nnca(nnc+j)=chain(j,i)
			lnca(nnc+j)=lbchain(j,i)
		end do 
		nnc=nnc+ichain(i)
	end if 
end do 

nchain=tnc ! corrected number of chain
write(*,*) "Number of Chains is: ", nchain 

do i=1,nchain
	write(*,*) "ichain: ",i,ichain(i) 
	write(*,*) "chain: ",i
	do j=1,ichain(i)
		write(*,"(xi0)",ADVANCE="NO") chain(j,i)
	end do 
	write(*,*) 
end do 
write(*,*) "complete read itp file"
!nchain,chain,ichain
END SUBROUTINE 

