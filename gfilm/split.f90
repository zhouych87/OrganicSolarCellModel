      program main 
!!ccccccc molecule-atom-distance-basis-core cccccccccccccccccc
      implicit none 
      integer::nam(1000),nm,i,j,k,l,m,n,count(1000),ii,rt,ri,tm,iin,iim  !number of atoms in molecule, and numbers of molecules
      real::c(3,1000,1000),dc(3,1000),cen(3,1000),dis(3),cutoff,dtc(3,1000) !c(im,ia,3)
      real::tdis,cell(3),kct,tdis1,tdis2(3)
      character(2)::elemt(1000,1000),delemt(2000) ! ,ALLOCATABLE
      character(20)::flname,mnum,mnum2,cnum,cnum2,flnm
      integer::st(1000,100),stt(1000,1000),km(1000),ncpu,tmp,tmp2 !,pair(1000,1000,2)
      character(20)::sets 
      character(100)::cmd
      integer::nhomo,mhomo
! 
!      write(*,*) 'Please input the coordinate file name'
!      read(*,*) flname
      open(100,file="in",status='old') 
      read(100,*) flname 
      read(100,*) ncpu 
      read(100,*) cutoff 
      read(100,*) 
      read(100,"(a20)") sets
      close(100)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     Read molecules coordinate           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dc=0
      call readgro(flname,cell,c,elemt,nm,nam,km) ! coordinate in nm
      write(*,"(3(xf12.3))") cell(1:3)
      Write(*,*) "Take the molecule separated by the boundary condition together" 
      do m=1,nm 
         do n=1,nam(m) 
           if (n/=1) then 
           do i=1,3
            if (abs(c(i,n-1,m)-c(i,n,m))>5) c(i,n,m)=c(i,n,m)+cell(i)*(c(i,n-1,m)-c(i,n,m))/abs(c(i,n-1,m)-c(i,n,m))
           end do !i
           end if 
        end do ! n
     end do ! m 
!##################################################
      !     calculate center coord of molecule    !
!##################################################
      kct=0;        cen=0;
      write(*,*) "molecules: ", nm 
      do m=1,nm 
 !        write(*,*) "atoms in molecule: ", m ,nam(m)
         do n=1,nam(m) 
             cen(:,m)=cen(:,m)+c(:,n,m)
         end do 
         cen(:,m)=cen(:,m)/nam(m)
!         Write(*,*) m , "central calculation finished" 
      end do 
!Write(*,*) "central calculation finished" 


!##################################################
!    select  dimer 
!##################################################
      dis=10000 
      open(131, file = "connection.dat", status = 'replace') 
!      write(131,*) nm ! number of molecules 
do m=1,nm 
    count(m)=0
    do i=-1,1  
    do j=-1,1  
    do k =-1,1  

       do n=1,nm 
         dis(1)=cen(1,n)+i*cell(1)-cen(1,m) 
         dis(2)=cen(2,n)+j*cell(2)-cen(2,m) 
         dis(3)=cen(3,n)+k*cell(3)-cen(3,m)
         tdis=dis(1)**2+dis(2)**2+dis(3)**2
       if (n/=m .and. tdis<cutoff**2 ) then 
!             write(*,*) "cell i j k",m,n,i,j,k,sqrt(tdis)
             tdis1=1000
             do iin=1,nam(n)
             do iim=1,nam(m)
                   tdis2(1)=c(1,iin,n)+i*cell(1)-c(1,iim,m) 
                   tdis2(2)=c(2,iin,n)+j*cell(2)-c(2,iim,m) 
                   tdis2(3)=c(3,iin,n)+k*cell(3)-c(3,iim,m)
                   tdis=tdis2(1)**2+tdis2(2)**2+tdis2(3)**2
                  if (tdis<tdis1) tdis1=tdis  ! find minimum to tdis1 
             enddo !iin
             end do !iim
!##################################################
!    selection based on nearest atom distance
!##################################################
            if (tdis1 < 0.250) then 
                count(m)=count(m)+1
                write(131,"(5(xi0),4(xf10.4))") m, i,j,k,n,(dis(l),l=1,3),sqrt(tdis1)             
            end if 
            
         end if 
         end do !! iin
         
    end do !k 
    end do !j
    end do !i
!    write(*,*) "Finished", m 
end do !m
    close(131)

!##################################################
!    Generate rw input file
!##################################################

      open(131, file = "connection.dat", status = 'old') 
      open(191, file ='grp', status = 'replace') !gather the results from integral calculation 
      write(191,"(a)") "#!/bin/bash "

      open(211, file ='tcon', status = 'replace') !  information to run random walk simulation
      write(211,'(i0,2xa,2xe15.5)') nm,'1000000 ',kct ! information to run random walk simulation 

      st=1 ! represent all the fragment is not calculated
      stt=1 ! dime not calculated 
!      pair=99 
      rt=10000 ! initial time > 32, to create a new run scripts
      ri=1

do m=1,nm 
!##################################################
!    Generate list of job need to run 
!##################################################
  if (rt>32) then 
      write(flname,"(a,i0)") "run",ri 
      open(181, file =flname, status = 'replace') !command to run g09
      write(181,"(a)") "#!/bin/bash "
      rt=0; 
  end if 
  
   write(211,'(i0,xi0,xi0,3(xe13.5))') m,count(m),km(m),(1.0E-9*cen(m,i),i=1,3) ! tcon 
   if (count(m)==0) goto 330 ! if there is no neighbouring molecules then stop write 
   
!##################################################
!    Fragment 1 
!##################################################
   i=0; j=0; k=0 
   call ncell(tmp,i,j,k)
   write(mnum,'(i0,a,i0,a,i0)') m,"m",tmp,"r",ri ! 
!   write(*,"(a,xxi0)") "st",st(m,tmp)
!   if (st(m,tmp)==1) then 
 !     write(*,*) "test" 
      write(flname,'(a,a)') trim(mnum),'.com' 
  !    write(*,'(a,a)')  'g09  ',trim(flname) !, ".log"  !181 alone  
      write(181,'(a,a)')  'g09  ',trim(flname) !, ".log"  !181 alone  
      call writegjf(flname,nam(m),c(:,:,m),elemt(:,m),ncpu,sets)
      call getmo(nam(m),elemt(:,m),mHOMO)
      rt=rt+1
!      st(m,tmp)=0 
!   end if
   write(mnum,'(2(a))') trim(mnum),".log" ! 
   
!##################################################
!    Fragment 2: neighbours 
!##################################################
   do l=1,count(m)
       read(131,*) tm,i,j,k,n !!!! cell 1 2 3 and molecules 
       if (tm /= m) goto 340 
       call ncell(tmp2,i,j,k); 

       ! if the dimer has already calculated then, not need to calculate again
    if (stt(m,n).eq.0) then
!          write(flnm,'(2(i0,a))') m,'t',n,'.log' 
!         if (pair(m,n,1)==0 .or. pair(m,n,2)==0) then 
!         write(*,*) "stt is not the same with pair","m n tmp:",stt(m,n),m,n,tmp, tmp2
!         end if 
!         write(mnum,'(i0,a,i0,a)') m,"m",pair(m,n,1),".log" ! ! 
!         write(mnum2,'(i0,a,i0,a)') n,"m",pair(m,n,2),".log" ! ! 
         write(flnm,'(2(i0,a))') m,'t',n,'.log' 
!         write(181,"(a,xa15,xa15,xa15,xi0,xi0,xa,a,i0)")  "calv", mnum,mnum2,flnm, km(m),km(n),">>", "v.out",ri 
!         write(191,"(3(a))") "grep '",flnm, "tmp.out >> v.out" 

    else if(stt(n,m).eq.0) then 
!        if (pair(n,m,1)==0 .or. pair(n,m,2)==0) then 
!        write(*,*) "stt is not the same with pair","m n tmp:",m,n,tmp, tmp2
!        end if 
         write(flnm,'(2(i0,a))') n,'t',m,'.log' 
!         write(mnum,'(i0,a,i0,a)') n,"m",pair(n,m,1),".log" !  ! 
!         write(mnum2,'(i0,a,i0,a)') m,"m",pair(n,m,2),".log" ! ! 
!        write(181,"(a,xa15,xa15,xa15,xi0,xi0,xa)") "calv", mnum,mnum2,flnm, km(n),km(m),">> v.out"
!          write(191,"(3(a))") "grep '",flnm, "tmp.out >> v.out"

    else
       
        !fragment 2 
        write(mnum2,'(i0,a,i0,a,i0)') n,"m",tmp2,"r",ri! 

!        if ((st(n,tmp2)==1) ) then ! all of them not calculated 
             write(flname,'(a,a)') trim(mnum2),'.com' 
             write(181,'(a,a)')  'g09  ',trim(flname)  !,".log"  !181 alone 
             dtc(1,1:nam(n))=c(1,1:nam(n),n)+i*cell(1) 
             dtc(2,1:nam(n))=c(2,1:nam(n),n)+j*cell(2) 
             dtc(3,1:nam(n))=c(3,1:nam(n),n)+k*cell(3) 
         
             call writegjf(flname,nam(n),dtc,elemt(:,n),ncpu,sets)
             call getmo(nam(n),elemt(:,n),nHOMO)
              rt=rt+1
!             st(n,tmp)=0
!         end if

!##################################################
!    Dimers
!##################################################
         write(mnum2,'(a,a)') trim(mnum2), ".log"! 
         write(flnm,'(2(i0,a))') m,'t',n,'.log' 
         write(flname,'(2(i0,a))') m,'t',n,'.com ' 
         write(181,'(a,a)')  'g09  ',flname !, ".log"  !181 dimer  

         dc=0; delemt="NO"
         dc(:,1:nam(m))=c(:,1:nam(m),m);           
         delemt(1:nam(m))=elemt(1:nam(m),m)
         
          dc(:,(1+nam(m)):(nam(m)+nam(n)))=dtc(:,1:nam(n))
         delemt((1+nam(m)):(nam(m)+nam(n)))=elemt(1:nam(n),n)

         call writegjf(flname,(nam(m)+nam(n)),dc,delemt,ncpu,sets) 
          rt=rt+4
!         write(181,"(a,xa15,xa15,xa15,xxi0,xxi0,xa)") "calv", mnum,mnum2,flnm, km(m),km(n),">> v.out"
         write(181,"(a,xa15,xa15,xa15,xi0,xi0,xa,a,i0)")  "calv", mnum,mnum2,flnm, int(mhomo/2),int(nhomo/2),">>", "v.out",ri 
!          write(181,"(a)")  "wait"
         write(181,"(3(xa))") "sleep 30; rm",mnum2,flnm
         stt(m,n)=0 
     end if 

     write(191,"(3(a))") "grep ' ",trim(flnm), " ' tmp.out >> v.out" 
   end do!n
 write(181,"(2(xa))") " rm",mnum 
 
330 if (rt>32) then 
      close(181)
     ri=ri+1 
   end if 
end do !m


close(131)
if (rt<=32) close(181)
close(191)
close(211) 
stop 

340 write(*,*) "tm is not equal to m", tm,m 
close(131)
if (rt<=32) close(181)
close(191)
close(211) 
stop
      end program ! main



subroutine ncell(n,i,j,k)
  implicit none
  integer::i,j,k,n
    n=(i+1)*3**2+(j+1)*3**1+(k+1)*3**0+1
  return 
end subroutine 
