# Expanded same-trajectory source of truth

Canonical version: `expanded_same_trajectory_v1`.

- Discovered raw trajectories: 81 across 16 proteins.
- Primary complete cohort: 73 across 16 proteins.
- Late retained: 17.
- Late non-retained: 56.

The master CSV retains every discovered trajectory. Use
`eligible_primary_cohort == True` for the primary 18--20 ns to 70--100 ns
same-trajectory analysis. Never infer eligibility from non-null labels alone;
use the explicit flag and audit `primary_exclusion_reason`.

Outcome: complete nominal 0--100 ns DCD coverage and median corrected ligand
pose RMSD below 3 A over (70,100] ns.

This dataset is retrospective. It is the row-level numerical source of truth,
not evidence of prospective validation or experimental binding.
