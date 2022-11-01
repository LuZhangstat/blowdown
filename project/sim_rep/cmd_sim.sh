cd ../..
for (( var=1; var<=2; var++ ))
do
	echo "run $var-th sim"
	R CMD BATCH -$var ./project/sim_rep/sim_rep.R ./results/sim/sim_rep$var
done

