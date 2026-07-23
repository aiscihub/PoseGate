"""Portable evaluation of frozen same-trajectory forecasting policies."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from .config import Policy
from .exceptions import ConfigurationError


@dataclass(frozen=True)
class PolicyDecision:
    forecast: str
    recommendation: str
    uncalibrated_retention_score: float
    decision_threshold: float
    score_calibrated: bool
    standardized_features: Mapping[str, float]
    linear_predictor: float


def _sigmoid(value: float) -> float:
    if value >= 0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def evaluate_policy(
    policy: Policy, measurements: Mapping[str, float]
) -> PolicyDecision:
    """Evaluate frozen scaler and logistic coefficients without joblib."""
    missing = [name for name in policy.feature_names if name not in measurements]
    if missing:
        raise ConfigurationError(f"policy measurements are missing: {missing}")

    standardized: dict[str, float] = {}
    linear_predictor = policy.intercept
    for name, mean, scale, coefficient in zip(
        policy.feature_names,
        policy.scaler_mean,
        policy.scaler_scale,
        policy.coefficients,
        strict=True,
    ):
        value = float(measurements[name])
        z_value = (value - mean) / scale
        standardized[name] = z_value
        linear_predictor += coefficient * z_value

    score = _sigmoid(linear_predictor)
    retained = score >= policy.decision_threshold
    return PolicyDecision(
        forecast="POSE_RETAINED" if retained else "POSE_NONRETAINED",
        recommendation="CONTINUE" if retained else "CANDIDATE_FOR_PAUSE",
        uncalibrated_retention_score=score,
        decision_threshold=policy.decision_threshold,
        score_calibrated=policy.score_calibrated,
        standardized_features=standardized,
        linear_predictor=linear_predictor,
    )
