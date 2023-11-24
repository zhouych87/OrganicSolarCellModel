# Organic Solar Cell Model
Fourier Transform on the Distance Matrix for the Fast Estimation of Electronic Couplings in Organic Semiconductors
#Due to the storage limit, large file such as cip.dat and other descriptors files are removed. Please generate descriptors before your use.
#read more on the descriptors please visit: https://pubs.acs.org/doi/10.1021/acs.jcim.3c00786
#pubs.acs.org/jcim Article
Intermolecular 3D-MoRSE Descriptors for Fast and Accurate Prediction of Electronic Couplings in Organic Semiconductors
J. Chem. Inf. Model. 2023, 63, 5089−5096
#The program to generate descriptor is https://github.com/zhouych87/OrganicSolarCellModel/tree/master/desfftmpi

#how to generate CIP data

desfftmpi tet.itp geo cn 0.1 tet.gro 100

#how to get standard 3D MoRSE descriptors

desfftmpi tet.itp 3dm snsr 1.0 tet.gro 32

#How to get get intermolecular 3D MoRSE descriptors

desfftmpi tet.itp d3d snsr 1.0 tet.gro 32

#how to get conjugated 3D MoRSE descriptors

desfftmpi tet.itp 3dpi snsr 1.0 tet.gro 32

#how to get conjugated intermolecular 3D MoRSE descriptors

desfftmpi tet.itp d3dpi snsr 1.0 tet.gro 32

#how to get FTDM descriptors

desfftmpi tet.itp d3dpi cn 0.1 tet.gro 36

#More descriptors please refer to https://github.com/zhouych87/OrganicSolarCellModel/tree/master/desfftmpi

#How to reproduce Fig3 

cd tet
python trainsize.py > trainsize.log
grep "CIP" trainsize.log > trainCIP.log
grep "UM" trainsize.log > trainUM.log
grep "UIM" trainsize.log > trainUIM.log
grep "FTDM" trainsize.log > trainFTDM.log
python plttsize.py 

size_tet_MAE.png will be generated. It is Fig3.

#How to reproduce Fig4 

python transfer1.py
python transfer2.py

#How to reproduce Fig5 

python filmfitplot.py

# Please let me know if there was any bug. Thank you! 
My email: zhouych87@gmail.com.
