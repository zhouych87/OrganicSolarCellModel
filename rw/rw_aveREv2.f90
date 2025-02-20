      program main 
      use MPI
      implicit none
!cccc mpif90 -limf -lm rw.f CCCCCC
      Real(kind=8),allocatable ::std(:),td(:,:),tt(:),p(:,:),conk(:,:),cen(:,:)
      real(kind=8)::ttt,kct,tre(4,2),re(5)
      Real(kind=8),allocatable ::dp(:,:),econk(:,:),hconk(:,:),ehconk(:,:),heconk(:,:),conc(:,:,:)
      Real(kind=8)::pp,pi,tstd,avstd,avttt,tttt,ttstd
      Real(kind=8)::rd,sij,hi,hj
      real::start,finish
      Integer:: eh,eunit
      character(5)::tmp 
      character(12)::na
      integer::i,j,jj,t,tp,mn,tht,ht,ti
      integer,allocatable ::st(:),tj(:),times(:),ttimes(:),tcon(:),con(:,:),km(:)
      integer :: myid,ierr,npcs,status(MPI_STATUS_SIZE),ns
    
      call cpu_time(start)
      call random_seed ()
      call MPI_INIT( ierr )     
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )     
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npcs, ierr )
!      call init_random_seed(myid)

      pi=3.1415926
!!!!! read inrw
      open(111,file='inrw',status='old')
      read(111,*,err=401)  tht,eh, ns, eunit ! tht record interval; electron or hole: 1-electron, 2-hole; simulation times; Unit of E, 1 is eV, 2 is hatree
      if (eh .ne. 1 .and. eh .ne. 2 ) goto 401
      if (eunit .ne. 1 .and. eunit .ne. 2 ) goto 401

      read(111,*,err=401)  (tre(i,1),i=1,4) ! energies of  kind one: 00 0+ ++ +0
      !no-charge struct without charge; nocharge-struct-charged; charged-struct-charged;charged-struct-without charge
      read(111,*,err=401)  (tre(i,2),i=1,4)  !energies of  kind two: 00 0+ ++ +0
      if (eh==2) then !! for hole 
          read(111,*,err=401)  (tre(i,1),i=1,4) 
          read(111,*,err=401)  (tre(i,2),i=1,4)   
       end if 
      if (eunit==2) tre=27.2*tre 
       
       goto 301
401   write(*,*) "Error read inrw, the correct format should be like this:"
      write(*,*) "Line 1: 10000 2 200 2 # record interval; 1-electron, 2-hole; simulation times; Unit of E, 1 is eV, 2 is hatree"
      write(*,*) "Line 2: -210022.5434 -210024.872 -210024.9504 -210022.4613 !!! drcn7t e 00 0- -- -0"
      write(*,*) "Line 3: -89347.67184 -89349.55312 -89349.71027 -89347.52224 !!!pcbm e 00 0- -- -0"
      write(*,*) "Line 4: -210022.5434 -210016.9371 -210017.0408 -210022.4354 !!! drcn7t h  00 0+ ++ +0"
      write(*,*) "Line 5: -89347.67184 -89341.08069 -89341.14996 -89347.60429 !!!pvbm h  00 0+ ++ +0"
      close(111)
      stop
301   close(111)
      
      

if (myid .eq. 0) then 
      if (eh==1) then 
           open(11,file='emblty.out',status='replace') 
           open(191,file='epath.out',status='replace') 
           open(100,file='e_occ_times.out',status='replace') 
      else if (eh==2) then 
           open(11,file='hmblty.out',status='replace') 
           open(191,file='hpath.out',status='replace')
           open(100,file='h_occ_times.out',status='replace')
      else 
      write(*,*) "eh must be 1(electron) or 2(hole), stop"
          stop
      end if 
end if 
      allocate(std(ns),td(ns,3),tt(ns),st(ns),tj(ns))
      
      if (tht<1000) then
        if (myid==0) write(*,*) 'the time is too short!!!!!!'
        ht=tht
      else 
         ht=tht
      end if
!!! reorganisation energy for self-exchange
      re(1)=tre(4,1)+tre(2,1)-tre(3,1)-tre(1,1)! 1 to 1
      re(2)=tre(4,2)+tre(2,2)-tre(3,2)-tre(1,2)! 2 to 2 
      re(5)=0.5*(re(1)+re(2))
!!! free energy change  
      re(3)=tre(1,1)+tre(3,2)-tre(3,1)-tre(1,2) !  1 to 2 final - inital
      re(4)=-re(3) ! 2 to 1 
      if (myid .eq. 0) then
         write(*,*) 'RE 1',re(1)
         write(*,*) 'RE 2',re(2)
         write(*,*) 'DE 1 to 2',re(3)
         write(*,*) 'DE 2 to 1',re(4)
    end if 
!cccccc read in  tcon    cccccccccccccc
      open(112,file='tcon',status='old')
      read(112,*) mn! number of molecules, time of simulation, prefactors
      
      allocate(times(mn),ttimes(mn),tcon(mn),con(mn,50),km(mn))
      allocate(p(mn,50),conk(mn,50),cen(mn,3))
      allocate(dp(mn,50),econk(mn,50),hconk(mn,50),ehconk(mn,50),heconk(mn,50),conc(mn,50,3))
      
      times=0;ttimes=0 
      dp=0;tcon=0;cen=0
      
      do i=1,mn
        read(112,*) ti,tcon(i),km(i),(cen(i,j),j=1,3) 
!        write(*,*) i,ti,tcon(i),km(i),(cen(i,j),j=1,3) 
         if (ti/=i) then 
            write(*,*) "molecule not pair in tcon", i,ti 
            stop 
         end if 
      end do 
      close(112)

!ccccccc read dimer detazyx kct tci and calc hopping rate cccccc
      open(120,file='v.out',status='old') 
      open(121,file='connection.dat',status='old')
      open(122,file='coupling.dat',status='replace')
      
      if (myid==0) write(122,*) "m1 m2  x y z dij  LUMO-LUMO HOMO-HOMO LUMO-HOMO  HOMO-LUMO"
!      do i=1,mn
!         write(*,*) i,tcon(i),km(i)
!      end do 
!      read(*,*)
      P=0

             
if (myid .eq. 0) write(*,*) "   t(s)               l(cm^2) "
if (myid .eq. 0) write(11,*) "  t(s)               l(cm^2)  "


do i=1,mn
   pp=0
!             write(*,"(2(xi0),f12.4)")  i, tcon(i) ,pp
 if (tcon(i).ge.1) then
         do j=1,tcon(i)

!             write(*,*) "Begin molecule",i,j,size(p)
             read(121,*) t,t,t,t,con(i,j),(conc(i,j,jj),jj=1,3)
!             write(*,"(3(xi0),xa,xi0,3(xf12.4),xi0)") t,t,t,"connect ",con(i,j),(conc(i,j,jj),jj=1,3)
!              read(*,*)
             read(120,*) tmp,tmp,tmp,tmp,tmp, na, econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
!             write(*,"(3(xa),4(xf12.4),xi0)") tmp, na, 'v.out',econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
            if (myid==0) write(122,"(2(xi0),8(xf15.5))") i,con(i,j),(conc(i,j,jj),jj=1,3),&
 & sqrt(conc(i,j,1)**2+conc(i,j,2)**2+conc(i,j,3)**2),econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
 !              write(*,"(a,3(xi0))")  "test1 ",i, tcon(i) 
            
            !conc(i,j,:)= conc(i,j,:)*(1.0e-9)  ! from nm to m  
                          
             if (eh==1) then 
                  conk(i,j)=econk(i,j)*0.001
             else if (eh==2) then 
                  conk(i,j)=hconk(i,j)*0.001
             end if 
           !  if (myid==0) write(333,'(a,4(2xf8.4))') 'sij hi hj Jij',sij,hi,hj,conk(i,j)
           !  conk(i,j)=(conk(i,j)-0.5*(hi+hj)*sij)/(1-sij*sij)  !v=(J-S(HI+HJ)/2)/(1-S**2)
           !  if (myid==0) write(333,'(2xf8.4)') conk(i,j)
            ! conk(i,j)=(6.58212E-16)**2/(2*kct*conk(i,j)**2) !!!H**2/2*V**2*kct !!!
            ! for Marcus theory ,kct is reorganisation energy
! conk(i,j)=sqrt(kct*0.026/pi) * exp((kct+dg)*(kct+dg)/(0.104*kct))* (6.58212E-16)/conk(i,j)**2  !1/r MARCUS THEORY! 
            if (km(i)== 1 .and. km(con(i,j))==1) then 
                conk(i,j)=sqrt(re(1)*0.026/pi)* exp(re(1)/0.104)/(conk(i,j)**2) ! reciprocal of marcus rate, h= * (6.58212E-16)
                
            else if (km(i)== 2 .and. km(con(i,j))==2) then 
                conk(i,j)=sqrt(re(2)*0.026/pi)* exp(re(2)/0.104)/(conk(i,j)**2) ! reciprocal of marcus rate
                
            else if (km(i)== 1 .and. km(con(i,j))==2) then 
            conk(i,j)=sqrt(re(5)*0.026/pi)* exp((re(5)+re(3))*(re(5)+re(3))/(0.104*re(5)))/(conk(i,j)**2) ! reciprocal of FGR rate

            else if (km(i)== 2 .and. km(con(i,j))==1) then 
            conk(i,j)=sqrt(re(5)*0.026/pi)* exp((re(5)+re(4))*(re(5)+re(4))/(0.104*re(5)))/(conk(i,j)**2) ! reciprocal of FGR rate
               
            else 
               write(*,*) "km is not right",i,j,km(i),km(con(i,j))
               stop 
            end if 
!            write(*,"(2(xf15.8))") conk(i,j),kct
            !write(*,"(a,xi0,4(xg12.4))") "j=",j ,pp, conk(i,j),1/conk(i,j),pp+1/conk(i,j)
           !  if (myid==0) write(*,110) con(i,j),(conc(i,j,jj),jj=1,3),conk(i,j)
             pp=pp+1/conk(i,j)
            !write(*,"(a,xi0,4(xg12.4))") "j=",j ,pp, conk(i,j),1/conk(i,j),pp+1/conk(i,j)
         end do
         ! write(*,"(a,3(xi0))")  "test2 ",i, tcon(i) 
         do j=1,tcon(i)
           dp(i,j)=1/(pp*conk(i,j))
         end do 
         ! write(*,*) 
         p(i,1)=dp(i,1)
!          write(*,"(a,4(xg12.4))") "p           1",p(i,1),conk(i,1),sqrt(kct*0.026/pi)*exp(kct/0.102)*(6.58212E-16)
        if (tcon(i).ge.2) then
         do j=2,tcon(i)
           p(i,j)=dp(i,j)+p(i,j-1)
!          write(*,"(i0,xa,xi0,2(xg12.4))") i,'tp ',j, p(i,j),conk(i,j) ! time of a hop
         end do 
         end if 
    end if 
!    write(*,"(a,3(xi0))")  "test3 ",i, tcon(i) 
end do
 
      close(120)
     close(121)
     close(122)

!     read(*,*) t 
     
!ccccccccc  init walk ccccccc 
!      open(10,file='rw.out',status='replace') 
!      if (myid==0) write(*,*) 'reading ok! init to random walk '


call init_random_seed(myid)
call random_number(rd)
      tt=0
      std=0
      st=int(rd*mn)
      tj=1
      
     do i=1,ns
      td(i,1:3)=cen(1,1:3)
      end do 

!      do i=1,ns
!19        call random_number(rd)
!        st(i)=ceiling(mn*rd)
!        if (tcon(st(i))==0) goto 19 
!       write(*,*) "st=",st(i)
!      end do 
!      read(*,*)
!do jj=1,5000
 
avttt=0;avstd=0
!do  while (avttt<=1.0e-6) !!! number of hopping,points to plot a figure!!!
do ti=1,500
       tstd=0;      ttt=0 ! before summary of j th simulation, let the initial value to be zero
       
         do 20 j=myid+1,ns,npcs !!! number of simulations !!!!
             do i=1,ht   !!!! save and calculate the distance !!!!
               !call init_random_seed(j)
                call random_number(rd)
                do tp=1,tcon(st(j))
                   if (rd<=p(st(j),tp)) then
                      tt(j)=conk(st(j),tp)+tt(j)  ! time of j th simulation 
                      td(j,:)=td(j,:)+conc(st(j),tp,:) ! displace of j th simulation 
                      st(j)=con(st(j),tp)
                      times(st(j))=times(st(j))+1 
                      exit
                   end if 
                end do !tp
!                  write(*,*) "hop 1 finished"
 !              if (tt(j)>=tj(j)*1.0e-9) then  
 !              write(*,*) "reach time",j, tt(j),tj(j)
 !                  tj(j)=tj(j)+1 
  !                 goto 21  
  !              end if 
             end do ! i ht 
             
21           std(j)=td(j,1)*td(j,1)+td(j,2)*td(j,2)+td(j,3)*td(j,3)
             tstd=tstd+std(j)   ! sum the distance running on one core
             ttt=ttt+tt(j)      ! sum the time running on one core 
20       continue  ! j 

         call MPI_ALLREDUCE(ttt,tttt,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(tstd,ttstd,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         
         call MPI_ALLREDUCE(times,ttimes,mn,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)
         
       if (myid .eq. 0) then
            avstd=ttstd/ns
            avttt=tttt/ns
!            write(10,130) avttt,avstd  (td(1,j),j=1,3),
            write(*,130)  avttt*(6.58212E-16),avstd*1.0e-14 ! nm^2 to cm^2  
            write(11,130)  avttt*(6.58212E-16),avstd*1.0e-14 ! nm^2 to cm^2  
            write(191,"(4(xg16.8),xi0)")  (td(1,j),j=1,3),td(1,1)**2+td(1,2)**2+td(1,3)**2,st(1)
     endif
         
      if (tt(j)>=1.0e+10) exit 
end do ! while t1 
if (myid .eq. 0 .and. eh==1) write(*,*) 'Electron mobility is: ',avstd*1.0e-14/(6*avttt*0.026*(6.58212E-16)),'cm^2/Vs'
if (myid .eq. 0 .and. eh==2) write(*,*) 'Hole mobility is: ',avstd*1.0e-14/(6*avttt*0.026*(6.58212E-16)),'cm^2/Vs'

       if (myid .eq. 0) then
      do j=1,mn
          write(100,"(2(2xg16.8))") ttimes(j)
      end do 
     end if 
     close(11)
     close(191)
     close(100)

deallocate(std,td,tt,st,tj)
      
      call MPI_Finalize(ierr)
!c      write(*,*)  ttt, tstd  
      call cpu_time(finish)
      print '("Time = ",f8.3," seconds.")',finish-start
130   format(2(2xg16.8))
110   format(2xi0,4(2xg12.5))
      END program main
  
      SUBROUTINE init_random_seed(myid)
            INTEGER :: i, n, clock,myid 
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
            CALL SYSTEM_CLOCK(COUNT=clock)
            seed =clock+37*(/(i-1,i=1,n)/)+1891202*myid+1
            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)
       END SUBROUTINE
