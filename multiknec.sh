


echo "hello knec-parallel~"

make clean
make

echo "finish building knec."

parameters_dir="parameter_list"
filelist=`ls $parameters_dir`


for file in $filelist
do  

    echo 'begin running knec for '$file
    outdir='Data_'$file
    mkdir $outdir

    screenfile='screen_output_'$file'.txt'
    :> $screenfile

    ./snec 'parameter_list/'$file > $screenfile && mv $screenfile $outdir & 
  

done








