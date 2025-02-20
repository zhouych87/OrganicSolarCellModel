      SUBROUTINE readgro(flname,cell,c,elemt,nm,nam,km)
      implicit none
      real::c(3,1000,1000),tc(3),cell(3)
      character(2)::elemt(1000,1000) ! ,ALLOCATABLE
      character(10)::flname,mname,tmname
      character(5)::telemnt
      integer::i,j,km(1000),nm,na,ia,l1,l2,tmp,nam(1000)

      open(111,file=flname,status='old')
      km=1;
      nm=1;
      ia=0;
      nam=0 
      elemt="XO" 
      c=0
      
      read(111,*)
      read(111,*) na  ! number of atoms
      do i=1,na
          read(111,"(a,a,i5,3(f8.3))") tmname,telemnt,tmp,(tc(j),j=1,3) 
!          write(*,*) tmname,telemnt,tmp,(tc(j),j=1,3) 
 !         if (i==1) tcen(km,im,:)=tc(:)
          if (i>1 .and. trim(tmname) /= trim(mname)) then  !mname is the last one 
             nm=nm+1; ! sequance of molecule 
             ia=0
 !            tcen(km,im,:)=tc(:)
             L1=len(trim(tmname))
             L2=len(trim(mname))
             if(tmname(l1-2:l1) /= mname(l2-2:l2)) then 
                 km(nm)=km(nm-1)+1  ! the kind of molecule 
            else 
                 km(nm)=km(nm-1)
             end if 
!             write(*,*) tmname,mname,tmname(l1-2:l1),mname(l1-2:l1), km(nm) 
          end if 
  !          tcen(km,im,:)=tcen(km,im,:)+tc(:)
          nam(nm)=nam(nm)+1   ! the element number of this kind molecule  
          ia=1+ia
          c(:,ia,nm)=tc(:)
          telemnt=trim(adjustl(telemnt))
          
          if (len(telemnt) == 1) then
              elemt(ia,nm)=telemnt 
          else 
          if (telemnt(2:2)=="1".or. telemnt(2:2)=="2" .or. telemnt(2:2)=="3" .or. telemnt(2:2)=="4" .or. telemnt(2:2)=="5" .or. & 
        &  telemnt(2:2)=="6" .or.  telemnt(2:2)=="7" .or. telemnt(2:2)=="8" .or. telemnt(2:2)=="9" ) then
             elemt(ia,nm)=telemnt(1:1) 
          else 
             elemt(ia,nm)=telemnt(1:2)
          end if 
          end if 
!         if(km(nm)==2) write(*,"(i8,xi0,xa,xi5,3(xf8.3))") km(nm),nm,elemt(nm,ia),ia,(c(nm,ia,j),j=1,3)
          mname=tmname
      end do 
      read(111,*) (cell(i),i=1,3)
      close(111)
      write(*,*) 'read in finished'
      return 
      END SUBROUTINE  

!!!g09 gjf
      SUBROUTINE writegjf(flname,an,coord,elemt,ncpu,sets)
      implicit none
      real::coord(3,1000)
      character(2)::elemt(1000)
      character(20)::flname
      character(20)::sets
      integer::i,j,an,ncpu
      
      open(320, file=flname,status = 'replace')
       write(320,'(a)') '%mem=8GB'
!       write(320,'(a)') '%nproc=1'
       write(320,'(a,i0)') '%nprocs=',ncpu 
       write(320,'(a)')   '%NoSave'
       write(320,'(3(a))') '# ',sets ,' iop(3/33=1,5/33=3) nosymm'
       write(320,*)
       write(320,'(a)') 'mol'
       write(320,*)
       write(320,'(a)') '0 1'
      do i=1,an
        write(320,'(a1,3(3xf10.4))') trim(adjustl(elemt(i))),(10*coord(j,i),j=1,3) ! from nm to A
      end do
      write(320,*) ' '
      write(320,*) ' '
      close(320)
      END SUBROUTINE 


!!!! get HOMO and LUMO
SUBROUTINE getmo(an,elemt,nHOMO)
integer::i,j,an,nhomo
character(2)::elemt(1000)
nhomo=0
do i =1,an
       if (trim(adjustl(elemt(i)))=="H") THEN ;  nhomo=nhomo+1 
       ELSE if (trim(adjustl(elemt(i)))=="C")THEN;nhomo=nhomo+6
       ELSE if (trim(adjustl(elemt(i)))=="N")THEN;nhomo=nhomo+7
       ELSE if (trim(adjustl(elemt(i)))=="O")THEN;nhomo=nhomo+8
       ELSE if (trim(adjustl(elemt(i)))=="F")THEN;nhomo=nhomo+9
       ELSE if (trim(adjustl(elemt(i)))=="Si")THEN;nhomo=nhomo+14
       ELSE if (trim(adjustl(elemt(i)))=="P")THEN;nhomo=nhomo+15
       ELSE if (trim(adjustl(elemt(i)))=="S")THEN;nhomo=nhomo+16
       ELSE if (trim(adjustl(elemt(i)))=="Cl")THEN;nhomo=nhomo+17
       ELSE if (trim(adjustl(elemt(i)))=="Br")THEN;nhomo=nhomo+35
       ELSE 
          write(*,*) "unkown element, stop",i,trim(adjustl(elemt(i))),nhomo,elemt(i)
          stop
      END IF 
end do

END SUBROUTINE 

