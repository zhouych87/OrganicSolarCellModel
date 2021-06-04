# Split the PBC thin film

# Compile. gfortran is required
make 

This is output a splitting.exe


# How to use it

1, Revise "in" file

tmp.gro               ! name of thin film grop file
12                    ! number of cores to run gaussian
1.5 0.5               ! cutoff distance (nm): center, nearest atoms
428 261               ! number of donor and accepters 
b3lyp/6-31g           ! basis, please check the generated file


2, Perform the compiled splitting.exe. This step will generate lots of "com", run* files.
"com" files are the gaussain input file.
run* are the scripts to run gaussian, which has seperate all gaussian into groups. You can run them one by one, or run them at once in several jobs. 


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
