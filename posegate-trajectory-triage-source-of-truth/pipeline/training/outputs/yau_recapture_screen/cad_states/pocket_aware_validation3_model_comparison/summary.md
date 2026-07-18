# Validation 3: pocket-aware state vs. simpler baselines

n=59 complexes, 15 proteins. y_stable100 base rate = 0.305.

                          n  balanced_accuracy  sensitivity  specificity  candidate_preservation  follow_up_fraction
approach                                                                                                            
existing_binary_label  59.0           0.464770     0.222222     0.707317                0.222222            0.271186
pose_only              59.0           0.568428     0.722222     0.414634                0.722222            0.627119
retention_only         59.0           0.537263     0.611111     0.463415                0.611111            0.559322
pocket_aware_combined  59.0           0.540650     0.666667     0.414634                0.666667            0.610169
