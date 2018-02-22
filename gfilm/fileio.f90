      SUBROUTINE readgro(flname,cell,c,elemt,nm,nam,km)
      implicit none
      real::c(1000,1000,3),tc(3),cell(3)
      character(7)::elemt(1000,1000) ! ,ALLOCATABLE
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
          c(nm,ia,:)=tc(:)
          elemt(nm,ia)=telemnt 
         
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
      real::coord(1000,3)
      character(7)::elemt(1000)
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
        write(320,'(a1,3(3xf10.4))') trim(adjustl(elemt(i))),(10*coord(i,j),j=1,3) ! from nm to A
      end do
      write(320,*) ' '
      write(320,*) ' '
      close(320)
      END SUBROUTINE 
