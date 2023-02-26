# split the PBC thin film

#compile gfortran is required
make 

This is output a split, copy the generated split into PATH folders. 


# How to use it
1, create a  in file


pen.gro    G54A7FF.itp          ! name of thin film grop file
24                    !core numbers
0.5               ! cutoff distance nm: nearest atoms

wB97XD/Def2TZVPP           ! method and basis set



2, Perform the compiled split. This step will generate lots of "com", run* files.
"com" files are the gaussain input file.
"run*" are the scripts to run gaussian, which has seperate all gaussian into groups. You can run them one by one, or run them at once in several jobs. 
"*.dat" are calculated discriptors

3, After splitting the film, you may use  the following command to find Gaussian calculation that not necessary. 
for i in {1..xxx}; do echo $i; grep -A1 g16 run$i |grep rm ; done

outputs:
1
2
3
....
186
187
  rm 442m14r187.log      
188
  rm 449m14r188.log      
189
  rm 451m14r189.log      

That means in the run187 document, the calculation of 442m14r187.com is not necessary. You can remove this line and the line before

g09 442m14r187.com
 rm 442m14r187.log
 
 4, perform Gaussian calculations. Please make sure g09 or g16 is in the PATH.
 for i in {1..100}
 do
 ./run$i
 done
