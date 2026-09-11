"""Portable evaluation of frozen same-trajectory forecasting policies."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from .config import Policy
from .exceptions import ConfigurationError


@dataclass(frozen=True)
class PolicyDecision:
    """One frozen rule applied to one measured prefix.

    Score fields are populated only by ``standardized_logistic`` policies.
    Margin fields are populated only by one-sided ``threshold_rule`` policies,
    which deliberately emit no probability of any kind.
    """

    forecast: str
    recommendation: str
    policy_type: str
    feature_values: Mapping[str, float]
    uncalibrated_retention_score: float | None = None
    decision_threshold: float | None = None
    score_calibrated: bool = False
    standardized_features: Mapping[str, float] | None = None
    linear_predictor: float | None = None
    stop_threshold_angstrom: float | None = None
    signed_stop_margin_angstrom: float | None = None


def _sigmoid(value: float) -> float:
    if value >= 0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def _feature_values(
    policy: Policy, measurements: Mapping[str, float]
) -> dict[str, float]:
    missing = [name for name in policy.feature_names if name not in measurements]
    if missing:
        raise ConfigurationError(f"policy measurements are missing: {missing}")
    values = {name: float(measurements[name]) for name in policy.feature_names}
    if not all(math.isfinite(value) for value in values.values()):
        raise ConfigurationError("policy measurements must be finite")
    return values


def _evaluate_standardized_logistic(
    policy: Policy, values: Mapping[str, float]
) -> PolicyDecision:
    """Evaluate frozen scaler and logistic coefficients without joblib."""
    standardized: dict[str, float] = {}
    linear_predictor = float(policy.intercept or 0.0)
    for name, mean, scale, coefficient in zip(
        policy.feature_names,
        policy.scaler_mean or (),
        policy.scaler_scale or (),
        policy.coefficients or (),
        strict=True,
    ):
        z_value = (values[name] - mean) / scale
        standardized[name] = z_value
        linear_predictor += coefficient * z_value

    score = _sigmoid(linear_predictor)
    threshold = float(policy.decision_threshold or 0.0)
    retained = score >= threshold
    return PolicyDecision(
        forecast="POSE_RETAINED" if retained else "POSE_NONRETAINED",
        recommendation="CONTINUE" if retained else "CANDIDATE_FOR_PAUSE",
        policy_type=policy.type,
        feature_values=dict(values),
        uncalibrated_retention_score=score,
        decision_threshold=threshold,
        score_calibrated=policy.score_calibrated,
        standardized_features=standardized,
        linear_predictor=linear_predictor,
    )


def _evaluate_threshold_rule(
    policy: Policy, values: Mapping[str, float]
) -> PolicyDecision:
    """Evaluate a one-sided early-stop screen.

    The rule fires only in the stop direction. Not firing is *not* a prediction
    that the pose is retained, so the forecast is ``DEFER`` rather than
    ``POSE_RETAINED``: the stop condition simply was not met. A value exactly
    equal to the threshold does not stop.
    """
    feature = policy.feature_names[0]
    value = values[feature]
    threshold = float(policy.stop_threshold_angstrom or 0.0)
    margin = value - threshold
    stop = margin > 0.0
    return PolicyDecision(
        forecast="POSE_NONRETAINED" if stop else "DEFER",
        recommendation="CANDIDATE_FOR_PAUSE" if stop else "CONTINUE",
        policy_type=policy.type,
        feature_values=dict(values),
        score_calibrated=False,
        stop_threshold_angstrom=threshold,
        signed_stop_margin_angstrom=margin,
    )


def evaluate_policy(
    policy: Policy, measurements: Mapping[str, float]
) -> PolicyDecision:
    values = _feature_values(policy, measurements)
    if policy.type == "standardized_logistic":
        return _evaluate_standardized_logistic(policy, values)
    if policy.type == "threshold_rule":
        return _evaluate_threshold_rule(policy, values)
    raise ConfigurationError(f"unsupported policy.type: {policy.type}")
