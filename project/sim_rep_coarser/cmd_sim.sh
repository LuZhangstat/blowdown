cd ../..
for (( var=1; var<=101; var++ ))
do
	echo "run $var-th sim"
	R CMD BATCH -$var ./project/sim_rep_coarser/sim_rep2.R ./results/sim2/sim_rep$var
done

