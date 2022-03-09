# We have two versions:
* rw_aveRE_CT. record time and distnace at certain period
* rw_aveRE     record time and distnace at centain steps

# The input file is inrw
```
100000 2 20000  # steps or periods to record a data; 1 for electron,2 for hole; 20000 number of simulations, 
-162874.6546 -162876.9996 -162877.1079 -162874.5371 !!! drcn5t e 00 0- -- -0
-89347.67184 -89349.55312 -89349.71027 -89347.52224 !!!pcbm e 00 0- -- -0
-162874.6546 -162868.6874 -162868.7941 -162874.5449  !!! drcn5t h 00 0+ ++ +0
-89347.67184 -89341.08069 -89341.14996 -89347.60429 !!!pvbm h 00 0+ ++ +0 
```
00 0- -- -0

Four set symbols, the first character of symbol is the charge of optimized the structure. The second is the charge of energy calculation

* 00 using the zero charge optimized structure to calculate the zero charge state energy
* 0- using the zero charge optimized structure to calculate the -1 charge state energy
* -- using the -1 charge optimized structure to calculate the -1 charge state energy
* -0 using the -1 charge optimized structure to calculate the zero charge state energy

# When the output line is not linear, most of the case are due to the error of reading coupling
Please check the v.dat and coupling.dat. v.out is the original coupling dat, and coupling.dat is the program readed dat please compare them to find the difference


# You may use below cmd to compare
cp v.out bv.out
sed -i "s/ Tll Thh Thl Tlh for  / /g" bv.out
sed -i "s/.log   / /g" bv.out 
sed -i "s/t/   /g" bv.out

The bv.out gives:  
```
m1 m2 vll vhh        vhl           vlh
 1   2    1.03489      2.24003      1.24722     -2.95822
 1   234    0.36956      0.57131     -0.29864     -1.29709
 1   291   -1.09224      1.47262      0.91710     -0.73112
 1   2    1.03489      2.24003      1.24722     -2.95822
 2   6   -4.27009      3.93086      7.91042     -1.90527
......
```

The coupling give 
```
 m1 m2         x               y               z               dij           vLUMO-LUMO     vHOMO-HOMO     vLUMO-HOMO  v   HOMO-LUMO
 1 2         0.79180        -0.16870         0.33030         0.87436         1.03489         2.24003         1.24722        -2.95822
 1 234         0.78080         0.02660        -0.63510         1.00683         0.36956         0.57131        -0.29864        -1.29709
 1 291        -0.87250         0.39360         0.02630         0.95753        -1.09224         1.47262         0.91710        -0.73112
 2 1        -0.79180         0.16870        -0.33030         0.87436         1.03489         2.24003         1.24722        -2.95822
```

Please make sure at the same line, m1 and m1 are the same in bv.out and coupling.out.

# If bv.out has more lines, please delete the dulipcate lines in v.out.

# If coupling has more lines, it v.out lose the data. It is caused by the failure in the QM calculations. You can calculate the missing coupling data bu hand.  

For example:
in coupling.dat 
```
286	115	0.4351	-0.3522	-1.1479	1.27712	-15.14306	-5.72938	-33.36647	-0.0554
286	164	0.335	1.2145	-0.3829	1.31676	0.94735	0.12212	2.58079	-0.39559
286	272	0.2736	0.0341	-0.8185	0.86369	-1.59779	0.50686	0.22062	0.2256
```

in bv.out
```
115	286	-0.25037	0.08672	0.04926	-0.61619
272	286	1.38739	-2.44986	-2.17782	-1.20501
```

It means the coupling of 164-286 dimer is missing.

```
zhouych@zhou0:~/rm$ grep "164t286.log" grp 
grep ' 164t286.log ' tmp.out >> v.out
grep ' 164t286.log ' tmp.out >> v.out
zhouych@zhou0:~/rw$ grep ' 164t286.log ' tmp.out
zhouych@zhou0:~/rw$ grep ' 164t286.log ' run*
run67:orbital=`g09log 164t286.log 1  | tail -1`
run67:echo 164t286.log $orbital  >> level.txt
run67: calv 164m14r67.log   286m14r67.log   164t286.log     383 647 >>v.out67
run67: sleep 10; rm 286m14r67.log        164t286.log        
```
Then you can calculate the 164-286 dimer coupling by hand:
The coupling is calculated by "calv 164m14r67.log   286m14r67.log   164t286.log     383 647 >>v.out67"
But before that, we need to calculate the 164m14r67.log   286m14r67.log   164t286.log

```
zhouych@zhou0:~/rw$ g09 164m14r67.com
zhouych@zhou0:~/rw$ g09 286m14r67.com
zhouych@zhou0:~/rw$ g09 164t286.com
zhouych@zhou0:~/rw$ calv 164m14r67.log   286m14r67.log   164t286.log     383 647 
 Tll Thh Thl Tlh for  1t234.log      0.36956      0.57131     -0.29864     -1.29709
```

add this line to the v.out at correspdong line(the line number is the same to the coupling data).Please check the m1 and m2 of above line and the below line in v.out and coupling.data is the same.



