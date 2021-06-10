! revised for the DOS calculation of OSC
! will output all levels of orbitals in the range of -7.0 eV -2.0 eV 
program g09log
  implicit none
  character(20):: fname
  character(27)::tmp1,mem,ncpu,chk
  character(79)::tmp,setting
  character(9)::val
  character(1)::arg2
  character(2)::elemnt(1000)
  integer::funct, nmo,lngth,fnmo,nline,i,j,nm,tmpi,nhm,nlm,hi,li 
  real:: hm(100),lm(100),crd(3,1000)
  logical:: omo,debug ! 0 false occupation 
  
  debug=.false. !
!  debug=.true.
  if (iargc()/=2) then 
     write(*,*) "input parameter missing,stop"!    ch1(ii)=c(ii,imol1)
     write(*,*) "please input the name of the gaussian output file, functions. "
     write(*,*) " FOR example: pcbm.log  1 " ! cl1(ii)=c(ii,imol1+1) !! add by Yecheng Zhou
     write(*,*) "0: obtain converge energy;"
     write(*,*) "1: obtain homo and lumo and its orbital number;"
     write(*,*) "2: write gjf from log, have to check the functional line"
     stop
  end if 
  
  call getarg(1,fname)
  call getarg(2,arg2)
  read(arg2,*) funct ! 0 convg; 1 homo,lumo;3 log2gjf 

if (funct==0) then 
!   setting="grep 'Predicted change in Energy' "//fname
!   fname2="|sed 's/\=/\ \ /'| sed 's/D/E/'"
!   write(*,'(a,a)') trim(setting),fname2
   setting="grep 'Predicted change in Energy' "//trim(fname)//"|sed 's/\=/\ \ /'| sed 's/D/E/'"
   write(*,*) trim(setting)
  CALL SYSTEM(trim(setting))
  stop 
end if 

 open(unit=99,file=fname,status='old',err=9993)

if (funct==1) then ! obtain homo and lumo and its orbital number
  nmo=0 
 omo=.True.
 nline=0
 hm=0.d0;lm=0.0d0
 hi=1; li=1 
 do while(.true.)
    read(99,'(a79)',end=9994) tmp
    if(tmp(1:27) .eq.' Alpha  occ. eigenvalues --') then 
       lngth=len(trim(tmp))-26 
       nmo=int(lngth/10)
       do i=1,nmo 
          val=tmp(30+10*(i-1):30+10*i)
          read(val,*) hm(hi)
          if (debug) write(*,"(i0x,f10.6x)") hi,hm(hi)
          hm(hi)=hm(hi)*27.2114
          if (hm(hi) > -7.d0) hi=hi+1
       end do 
       nhm=hi-1

    else if ((omo .eqv. .True.) .and. (tmp(1:27).eq.' Alpha virt. eigenvalues --')) then 
       lngth=len(trim(tmp))-26 
       nmo=int(lngth/10)
       do i=1,nmo 
          val=tmp(30+10*(i-1):30+10*i)
          read(val,*) lm(li)
          if (debug) write(*,"(i0x,f10.6x)") li,lm(li)
          lm(li)=lm(li)*27.2114 
          if (lm(li) > -2.d0) then 
            nlm=li
            omo=.false.
          end if 
          li=li+1
       end do 
           
    else if(tmp(1:27) .eq. ' Normal termination of Gaus') then 
        exit
    end if 
 enddo


do i=1,nhm
    write(*,"(f10.6x)",advance="No") hm(i)
enddo

write(*,"(a)",advance="No") "gap"
do i=1,nlm
    write(*,"(f10.6x)",advance="No") lm(i)
enddo
write(*,"(a)") " "

  stop

else if (funct==2) then ! funct=2, write gjf from geometry relaxation  log 
   omo =.false.
  do while(.true.)
    read(99,'(a)',end=9994) tmp
    !write(*,*) trim(tmp) 
    
    if (omo .eqv. .false.) then 
        if (tmp(1:5) .eq. ' %mem') then 
          mem=tmp(2:20) 
        else if (tmp(1:6) .eq. ' %npro') then 
          ncpu= tmp(2:20) 
        else if (tmp(1:2) .eq. ' #') then 
          setting= tmp(2:len(tmp))
          chk=fname(1:(len(trim(fname))-4))
          write(*,*) 'chk',chk,fname(1:(len(trim(fname))-4)),fname
          write(*,*) "setting ok"
        else if (tmp(1:19) .eq. ' Symbolic Z-matrix:') then 
            write(*,*) "going to read element and coordination"
            read(99,'(a)',end=9994) tmp
             i=0
             do while(.true.)
                i=i+1
                read(99,'(a)',end=9994) tmp
                if (tmp(1:2) .eq. '  '  ) then 
                     nm=i-1
                     write(*,*) 'Setting and element read finished, nm',nm 
                     goto 9990
                else 
                    read(tmp(1:2),"(2a)",err=9991,end=9995) elemnt(i)
                end if 
              end do 
         
9991          write(*,*) tmp
              write(*,*) 'error in reading element'
              stop 
              
9990          omo=.true.
              write(*,*) "setting all ok"
              write(*,*) mem,chk,ncpu  
        end if ! crd 
        
     else if (tmp(1:45) .eq. & 
     & '                          Input orientation:' .or. &
     & tmp(1:47) .eq. & 
     & '                         Z-Matrix orientation:' .or. &
     & tmp(1:47) .eq. & 
     & '                         Standard orientation:') then 
          write(*,*) 'Update coordination +1 '
   
          do i=1,4
             read(99,'(a)')
          end do 
          do i=1,nm
             read(99,*,end=9996) j,j,j,(crd(j,i),j=1,3)
          end do 
     else if (tmp(1:29) .eq. ' Number     Number       Type') then  ! update 20190103
          write(*,*) 'Update coordination +1 '
          read(99,'(a)')
          do i=1,nm
             read(99,*,end=9996) j,j,j,(crd(j,i),j=1,3)
          end do   
     end if 
     
    if(tmp(1:27) .eq. ' Normal termination of Gaus') then 
       Write(*,*) 'Gaussian terminated normally'
        goto 9992  
     end if 
   end do 
   write(*,*) "Abnormal termination of Gaussian!!!!!!!"
   
9992   write(*,'(a,a,a)') 'Going to write ', trim(chk), 'mo.gjf' 
   call writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
   write(*,*) 'writing finished'
   stop 
end if 
    
9993 write(*,*) "not exist"
9994 write(*,*) "Calculation is not finished or not succeed"
9995 write(*,*) "element"
9996 write(*,*) "crd"
write(*,*) tmp

   write(*,'(a,a,a)') 'Going to write ', trim(chk), 'mo.gjf'
   call writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
   write(*,*) 'writing finished'
   stop 

end program g09log 



      SUBROUTINE writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
      implicit none
      real::crd(3,1000)
      character(2)::elemnt(1000)
      character(30):: fname
      character(27)::mem,ncpu,chk
      character(79)::setting  
      integer::i,j,nm
      
      
      write(fname,'(a,a)') trim(chk), 'mo.gjf' 
      open(320, file=trim(fname),status = 'replace')
      write(320,'(a,i0)') trim(ncpu)
      write(320,'(a)')  trim(mem)
      write(320,'(a,a,a)') '%chk=',trim(chk),'.chk'
      write(320,'(a)') trim(setting)
      write(320,'(a)') ' '
      write(320,'(a)') 'MO generated'
      write(320,'(a)')  ' '
      write(320,'(a)') '0 1' ! should check it
      
      do i=1,nm
        write(320,'(a,3(3xf12.6))') elemnt(i),(crd(j,i),j=1,3) 
      end do
      write(320,*) ' '
      write(320,*) ' '
      close(320)
      END SUBROUTINE 
