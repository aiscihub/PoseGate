# Validation 3: protein-clustered significance

Cohort: 59 complexes, 15 proteins.

pred_existing_binary_label: mean=0.465, 95% CI=[0.362, 0.590]
pred_pose_only: mean=0.571, 95% CI=[0.470, 0.679]
pred_retention_only: mean=0.542, 95% CI=[0.450, 0.644]
pred_pocket_aware_combined: mean=0.540, 95% CI=[0.454, 0.632]

pose_only - existing_binary_label balanced-accuracy difference: mean=0.106, 95% CI=[-0.039, 0.244]  (includes 0 -- not significant)
pocket_aware_combined - pose_only balanced-accuracy difference: mean=-0.031, 95% CI=[-0.125, 0.043]  (includes 0 -- not significant)
