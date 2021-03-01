python dae_core_test.py > ~/logs/dae_core_test.log & 
BACK_PID=$!
wait $BACK_PID

python dae_noise_test.py > ~/logs/dae_noise_test.log & 
BACK_PID=$!
wait $BACK_PID

python sae_l_test.py > ~/logs/sae_l_test.log & 
BACK_PID=$!
wait $BACK_PID

python cae_lam_test.py > ~/logs/cae_lam_test.log &