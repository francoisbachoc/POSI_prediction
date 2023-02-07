#loop over launching R scripts in parallel
# parameters : 1: smallest index of script run[index].R to run
#              2: largest index of script run[index].R to run

path_=equicorrelated_protected_intercept_p9_n20_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=equicorrelated_protected_intercept_p9_n100_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=exchangeable_protected_intercept_p9_n20_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=exchangeable_protected_intercept_p9_n100_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=watershed_protected_intercept_p9_n20_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=watershed_protected_intercept_p9_n100_SCAD_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..



path_=equicorrelated_protected_intercept_p9_n20_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=equicorrelated_protected_intercept_p9_n100_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=exchangeable_protected_intercept_p9_n20_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=exchangeable_protected_intercept_p9_n100_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=watershed_protected_intercept_p9_n20_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..

path_=watershed_protected_intercept_p9_n100_MCP_1_withRes
echo ----------------------------------------------------------------------------------
echo $path_
echo ----------------------------------------------------------------------------------
cd $path_
for i in `seq 1 8`
do
     echo $i
     R CMD BATCH run"$i".R  &     >>  out$i.txt   
done
cd ..
