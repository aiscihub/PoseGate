# Physics-guided AI audit

This audit compares models using only corrected physical measurements available
by 20 ns. Evaluation uses outer leave-one-protein-out predictions. Every tuned
model performs hyperparameter selection using only the outer training proteins
through an inner protein-grouped validation loop.

## Main result

The transparent corrected-RMSD score remains the strongest ranking model on the
73-trajectory, 16-protein expanded cohort (AUROC 0.825). None of the nonlinear
or temporal AI models improves on it:

- endpoint logistic regression: AUROC 0.818;
- temporal logistic regression: AUROC 0.797;
- temporal random forest: AUROC 0.788;
- endpoint spline logistic model: AUROC 0.785;
- temporal histogram gradient boosting: AUROC 0.712.

The paired temporal-minus-endpoint logistic AUROC difference is -0.021; a
1,000-draw protein-clustered bootstrap interval is approximately [-0.115,
0.067]. Thus the current data do not support a ranking improvement from the
added temporal AI features.

## Operational observation

At the fixed 0.5 score threshold, endpoint logistic regression stops 40 runs,
including 3 of 17 retained trajectories, and misses 19 of 56 non-retained
trajectories. Temporal logistic regression stops 42 runs with the same 3 false
stops and misses 17 non-retained trajectories. This is a small retrospective
policy improvement (two additional non-retained runs stopped), not evidence of
better AUROC, and requires prospective validation before promotion.

## Interpretation

The negative model-complexity result is scientifically useful: the signal is
predominantly a low-dimensional physical persistence signal at the decision
checkpoint rather than a hidden nonlinear pattern recovered by flexible AI.
The physical corrected-RMSD score should remain the primary claim and baseline.
AI can be retained as a secondary resource-allocation layer whose incremental
value is judged by false-stop-controlled compute savings.

## Files

- `model_comparison_metrics.csv`: outer-LOPO AUROCs and clustered intervals.
- `nested_lopo_predictions.csv`: row-level out-of-fold scores.
- `nested_model_selections.csv`: configuration selected inside each outer fold.
