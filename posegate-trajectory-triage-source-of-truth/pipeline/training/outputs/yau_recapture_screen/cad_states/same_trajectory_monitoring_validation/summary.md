# Same-trajectory monitoring validation (correct design for 'extend or stop' decision)

Cohort: 58 complexes, 15 proteins. y_stable100 base rate = 0.293.

## State distribution
s20_state_sametraj
pocket_lost_instability         28
retained_stable                 16
pocket_retained_perturbation     9
metric_discordant_uncertain      5

## 2x2 table (pocket-retained-perturbation vs. pocket-lost-instability)
pocket_retained_perturbation: n=9 (8 proteins), stable rate=0.222 (2/9)
pocket_lost_instability: n=28 (12 proteins), stable rate=0.036 (pred_combined/28)
Odds ratio=7.714, Fisher p=0.1405
Protein-clustered bootstrap: rate diff 95% CI=(np.float64(-0.06896551724137931), np.float64(0.4444444444444444)), OR 95% CI=(np.float64(0.7775424768355109), np.float64(38.81818181818182))

## Model comparison (LOPO, protein-clustered bootstrap)
pred_pose_only: mean balanced_accuracy=0.765, 95% CI=[0.593, 0.882]
pred_retention_only: mean balanced_accuracy=0.659, 95% CI=[0.562, 0.750]
pred_combined: mean balanced_accuracy=0.762, 95% CI=[0.639, 0.847]
