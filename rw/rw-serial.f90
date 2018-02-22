      program main 
      implicit none
!cccc mpif90 -limf -lm rw.f CCCCCC
      Real(16)::std(40000),td(40000,3),tt(40000),ttt,kct,tre(2,4),re(4)
      Real(16)::dp(1000,50),econk(1000,50),hconk(1000,50),ehconk(1000,50),heconk(1000,50),conc(1000,50,3)
      Real(16)::pp,p(1000,50),pi,tstd,avstd,avttt,tttt,ttstd,conk(1000,50)
      Real(16)::rd,sij,hi,hj,cen(1000,3)
      Integer:: eh 
      character(5)::tmp 
      character(12)::na
      integer::tcon(1000),con(1000,50),i,j,jj,t,st(40000),tp,mn,tht,ht,ti,km(1000)
      integer :: ns,times(1000),np

      call random_seed ()
!      call init_random_seed(myid)

      dp=0
      tcon=0
      pi=3.1415926
      cen=0 
      times=0 
      open(111,file='inrw',status='old')
      read(111,*)  tht,eh, ns ,np ! tht record interval,electron-hole,simulation times
      read(111,*)  (tre(1,i),i=1,4) ! energies of  kind one: 00 0+ ++ +0
      ! no-charge struct without charge; nocharge-struct-charged; charged-struct-charged;charged-struct-without charge
      read(111,*)  (tre(2,i),i=1,4)  !energies of  kind two: 00 0+ ++ +0
      if (eh==2) then !! for hole 
          read(111,*)  (tre(1,i),i=1,4) 
          read(111,*)  (tre(2,i),i=1,4)   
       end if 
       
      close(111)
      
      if (tht<1000) then
       write(*,*) 'the time is too short!!!!!!'
        ht=tht
      else 
         ht=tht
      end if
      
      re(1)=tre(1,4)+tre(1,2)-tre(1,3)-tre(1,1) ! 1 to 1
      re(2)=tre(2,4)+tre(2,2)-tre(2,3)-tre(2,1) ! 2 to 2
      re(3)=tre(1,4)+tre(2,2)-tre(1,3)-tre(2,1) !  1 to 2
      re(4)=tre(2,4)+tre(1,2)-tre(2,3)-tre(1,1) ! 2 to 1 
      
!cccccc read in  tcon    cccccccccccccc
      open(111,file='tcon',status='old')
      read(111,*) mn ! number of molecules, time of simulation, prefactors

      
      do i=1,mn
        read(111,*) ti,tcon(i),km(i),(cen(i,j),j=1,3) 
!        write(*,*) i,ti,tcon(i),km(i),(cen(i,j),j=1,3) 
         if (ti/=i) then 
            write(*,*) "molecule not pair in tcon", i,ti 
            stop 
         end if 
      end do 
      close(111)

!ccccccc read dimer detazyx kct tci and calc hopping rate cccccc
      open(191,file='v.out',status='old') 
      open(121,file='connection.dat',status='old')
      open(122,file='coupling.dat',status='replace')
      
!      do i=1,mn
!         write(*,*) i,tcon(i),km(i)
!      end do 
!      read(*,*)
      P=0
      
do i=1,mn
   pp=0
!             write(*,"(2(xi0),f12.4)")  i, tcon(i) ,pp
 if (tcon(i).ge.1) then
         do j=1,tcon(i)

!             write(*,*) "Begin molecule",i,j,size(p)
             read(121,*) t,t,t,t,con(i,j),(conc(i,j,jj),jj=1,3)
!             write(*,"(3(xi0),xa,xi0,3(xf12.4),xi0)") t,t,t,"connect ",con(i,j),(conc(i,j,jj),jj=1,3)
!              read(*,*)
             read(191,*) tmp,tmp,tmp,tmp,tmp, na, econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
!             write(*,"(3(xa),4(xf12.4),xi0)") tmp, na, 'v.out',econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
             write(122,"(2(xi0),8(xf15.5))") i,con(i,j),(conc(i,j,jj),jj=1,3),&
 & sqrt(conc(i,j,1)**2+conc(i,j,2)**2+conc(i,j,3)**2),econk(i,j),hconk(i,j),ehconk(i,j),heconk(i,j)
 !              write(*,"(a,3(xi0))")  "test1 ",i, tcon(i) 
            
            conc(i,j,:)= conc(i,j,:)*(1.0e-9)  ! from nm to m  
                          
             if (eh==1) then 
                  conk(i,j)=econk(i,j)*0.001 ! max(abs(econk(i,j)*0.001),0.000000003 )
             else if (eh==2) then 
                  conk(i,j)=hconk(i,j)*0.001
            else 
                 write(*,*) "eh must be 1(electron) or 2(hole), stop"
                 stop
             end if 
           !  if (myid==0) write(333,'(a,4(2xf8.4))') 'sij hi hj Jij',sij,hi,hj,conk(i,j)
           !  conk(i,j)=(conk(i,j)-0.5*(hi+hj)*sij)/(1-sij*sij)  !v=(J-S(HI+HJ)/2)/(1-S**2)
           !  if (myid==0) write(333,'(2xf8.4)') conk(i,j)
            ! conk(i,j)=(6.58212E-16)**2/(2*kct*conk(i,j)**2) !!!H**2/2*V**2*kct !!!
            ! for Marcus theory ,kct is reorganisation energy
            if (km(i)== 1 .and. km(con(i,j))==1) then 
                kct=re(1)
            else if (km(i)== 2 .and. km(con(i,j))==2) then 
                kct=re(2)
            else if (km(i)== 1 .and. km(con(i,j))==2) then 
               kct=re(3)
            else if (km(i)== 2 .and. km(con(i,j))==1) then 
               kct=re(4)
            else 
               write(*,*) "km is not right",i,j,km(i),km(con(i,j))
               stop 
            end if 
!            write(*,"(2(xf15.8))") conk(i,j),kct

             conk(i,j)=sqrt(kct*0.026/pi)*exp(kct/0.102)*(6.58212E-16)/(conk(i,j)**2) 
            ! write(*,"(a,xi0,4(xg12.4))") "j=",j ,pp, conk(i,j),1/conk(i,j),pp+1/conk(i,j)
           !  if (myid==0) write(*,110) con(i,j),(conc(i,j,jj),jj=1,3),conk(i,j)
             pp=pp+1/conk(i,j)
            ! write(*,"(a,xi0,4(xg12.4))") "j=",j ,pp, conk(i,j),1/conk(i,j),pp+1/conk(i,j)
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
!          write(*,"(i0,xa,xi0,2(xg12.4))") i,'tp ',j, p(i,j),conk(i,j)
         end do 
         end if 
    end if 
!    write(*,"(a,3(xi0))")  "test3 ",i, tcon(i) 
end do
!     write(*,*) "read finished"
      close(191)
     close(121)
     close(122)

!     read(*,*) t 
     
!ccccccccc  init walk ccccccc 
!      open(10,file='rw.out',status='replace') 
!      if (myid==0) write(*,*) 'reading ok! init to random walk '

      open(191,file='path.out',status='replace') 
      tt=0
      std=0
      st=1 
      td=0 
      

call init_random_seed()
      do i=1,ns
19        call random_number(rd)
        st(i)=ceiling(mn*rd)
        if (tcon(st(i))==0) goto 19 
 !      write(*,*) "st=",st(i)
!        td(j,1:3)=cen(j,1:3)
      end do 
!      read(*,*)
!do jj=1,5000

!do  while (avttt<=1.0e-6) !!! number of hopping,points to plot a figure!!!
do ti=1,np
       tstd=0;      ttt=0 ! before summary of j th simulation, let the initial value to be zero
       
      do  j=1,ns !!! number of simulations !!!!
!          do i=1,ht   !!!! save and calculate the distance !!!!
          do while (tt(j) <= ti*1.0e-8)
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
                
         end do ! i ht 
             
             std(j)=td(j,1)**2+td(j,2)**2+td(j,3)**2
             tstd=tstd+std(j)   ! sum the distance running on one core
             ttt=ttt+tt(j)      ! sum the time running on one core 
      end do  ! j 

            avstd=tstd/ns
            avttt=ttt/ns
!            write(10,130) avttt,avstd  (td(1,j),j=1,3),
            write(*,130)  avttt,avstd
            write(191,"(4(xg16.8),xi0)")  (td(1,j),j=1,3),td(1,1)**2+td(1,2)**2+td(1,3)**2,st(1)
         
      if (avttt>=1.0e-6) exit 
end do ! 5000 points


      do j=1,mn
          write(100,"(2(2xg16.8))") times(j)
      end do 
      
     close(191)

      

130   format(2(2xg16.8))
110   format(2xi0,4(2xg12.5))
      END program main
  
      SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
            CALL SYSTEM_CLOCK(COUNT=clock)
            seed =clock+37*(/(i-1,i=1,n)/)
            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)
       END SUBROUTINE
