# Validation 2: 20 ns pocket-aware state vs. separate 100 ns structural outcome

S100 rule: stable if pose_rmsd_late < 3.0 A AND com_displacement_late < 5.0 A (protein-CA-aligned, PBC-corrected, contact-free, median over last 30% of BASE_100NS).

## 2x2 table
                              stable_100ns  unstable_100ns
20ns_state                                                
pocket_retained_perturbation             6              11
pocket_lost_instability                  7              20

pocket_retained_perturbation: n=17 (10 proteins), stable-outcome rate = 0.353 (6/17)
pocket_lost_instability:      n=27 (12 proteins), stable-outcome rate = 0.259 (7/27)

Odds ratio (Fisher): 1.558
Fisher exact p-value: 0.5205
Protein-clustered bootstrap stable-rate difference 95% CI: [-0.194, 0.367]
Protein-clustered bootstrap odds ratio 95% CI: [0.436, 6.890]
