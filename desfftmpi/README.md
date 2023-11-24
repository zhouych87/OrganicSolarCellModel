# This program wrote in fortran to generate descriptors
# intel mpi is required.
# how to compile

copy all file into a folder and then input make in the terminal. "desfftmpi" will be generated

# copy it to your PATH folders. make sure you can use is just type desfftmpi in terminal.

# how to use:

mpirun -n 24 desfftmpi tet.itp d3dpi cn 0.1 tet.gro 100

tet.itp: is the force field file to generate films. It is required to find conjugated atoms
tet.gro: is the atom cooridnation of film, which the format of GROMACS structure.
d3dpi: tells the program what descriptors will be generated. See below for details.
cn: is the core function form. See below for details.
0.1: is scale factor
100: It the dimension of descriptors, which is the range of S in 3D MoRSE.


##what descriptors will be generated:

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
        
        

        
#Core functions

        ttt=sl*s*r
        if (trim(method)=="snsr") then  ! standard 3D MoRSE
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
        else if (trim(method)=="cn") then   ! FTDM
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
