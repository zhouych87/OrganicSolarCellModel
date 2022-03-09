# Split the thin film PBC unit cell

## Compile 
gfortran is required. 
`make`

This will output a splitting.exe. 


## How to use it

### 1, Revise "in" file, which is the required input file

    tmp.gro               ! name of thin film grop file     
    12                    ! number of cores to run gaussian     
    1.5                   ! center cutoff distance (nm): the nearest atom distance cutoff is 0.5 nm by default.       
    b3lyp/6-31g           ! basis, please check the generated file    


### 2, Perform the compiled splitting.exe. This step will generate lots of "com", run* files.    
"com" files are the gaussain input file.    
run* are the scripts to run gaussian, which has seperate all gaussian into groups. You can run them one by one, or run them at once in several jobs. 


### 3, After splitting the film, you may use  the following command to remove Gaussian calculation that not necessary.     
    
    for i in {1..300}     
    do     
    line= \`grep -nA3 "g09 " run$i| grep rm | awk '{print $1}'| sed  's/-/ /g'\`      
    for j in $line      
    do     
    sline=\`echo " $j - 3 "| bc\`      
    sed -i "$sline,${j}d" run$i    
    done  
    done   
