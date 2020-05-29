

#for i in MTX/mult_dcop_03.mtx MTX/mult_dcop_01.mtx MTX/mult_dcop_02.mtx MTX/lp_fit2d.mtx MTX/bloweya.mtx MTX/lp_osa_07.mtx MTX/ex19.mtx MTX/brainpc2.mtx MTX/shermanACb.mtx MTX/cvxqp3.mtx MTX/case9.mtx MTX/TSOPF_FS_b9_c6.mtx MTX/OPF_6000.mtx MTX/OPF_3754.mtx MTX/c-47.mtx MTX/mhd4800a.mtx MTX/gen4.mtx MTX/Maragal_6.mtx MTX/aft01.mtx MTX/TSOPF_RS_b39_c7.mtx
#do
#    python3 ent.py -f $i -d 0  
#done


FILES="MTX/mult_dcop_03.mtx MTX/mult_dcop_01.mtx MTX/mult_dcop_02.mtx MTX/lp_fit2d.mtx MTX/bloweya.mtx MTX/lp_osa_07.mtx MTX/ex19.mtx MTX/brainpc2.mtx MTX/shermanACb.mtx MTX/cvxqp3.mtx MTX/case9.mtx MTX/TSOPF_FS_b9_c6.mtx MTX/OPF_6000.mtx MTX/OPF_3754.mtx MTX/c-47.mtx MTX/mhd4800a.mtx MTX/gen4.mtx MTX/Maragal_6.mtx MTX/aft01.mtx MTX/TSOPF_RS_b39_c7.mtx"

FILES1="MTX/mult_dcop_02.mtx MTX/lp_fit2d.mtx MTX/bloweya.mtx MTX/lp_osa_07.mtx MTX/ex19.mtx MTX/brainpc2.mtx MTX/shermanACb.mtx MTX/cvxqp3.mtx MTX/TSOPF_RS_b39_c7.mtx"

FILES="MTX/case9.mtx MTX/TSOPF_FS_b9_c6.mtx MTX/OPF_6000.mtx MTX/OPF_3754.mtx MTX/c-47.mtx MTX/mhd4800a.mtx MTX/gen4.mtx MTX/Maragal_6.mtx MTX/aft01.mtx"




#python3 ent.py -f "MTX/mult_dcop_03.mtx MTX/mult_dcop_02.mtx MTX/lp_fit2d.mtx MTX/bloweya.mtx MTX/lp_osa_07.mtx MTX/#ex19.mtx MTX/brainpc2.mtx MTX/shermanACb.mtx MTX/cvxqp3.mtx MTX/TSOPF_RS_b39_c7.mtx" \
#	-gpu -d 1 -t 32 -s fiji1.bin -se --temp SHUFFLED1 > out11 2> err1 & 

python3 ent.py -f "MTX/mult_dcop_03.mtx MTX/mult_dcop_01.mtx MTX/mult_dcop_02.mtx MTX/lp_fit2d.mtx MTX/bloweya.mtx MTX/lp_osa_07.mtx MTX/ex19.mtx MTX/brainpc2.mtx MTX/shermanACb.mtx MTX/cvxqp3.mtx MTX/case9.mtx MTX/TSOPF_FS_b9_c6.mtx MTX/OPF_6000.mtx MTX/OPF_3754.mtx MTX/c-47.mtx MTX/mhd4800a.mtx MTX/gen4.mtx MTX/Maragal_6.mtx MTX/aft01.mtx MTX/TSOPF_RS_b39_c7.mtx" \
	-gpu -d 0 -t 32 -s fiji.bin -se --temp SHUFFLED > out 2> err 
