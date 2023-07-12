dir="/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2"
cd $dir
for sample in {48..51}
do
if test -e controlfile.txt
then
rm controlfile.txt
dev_pbs_index=$sample
echo $sample
echo "$dev_pbs_index" > /media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/index.txt
cd $dir
./run_relernn_nhm.sh
echo "Sample $sample started"
else
sleep 2h
fi
done

#resultfiles=()
#'find . -type f -name "*BSCORRECTED.txt" -print | grep -i '\.txt$''
#for bscor_file in resultfiles
#do
# Rscript /media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/PlotReLernnResults.R $bscorfile
#done

