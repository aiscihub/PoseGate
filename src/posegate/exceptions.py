"""Domain-specific failures with actionable command-line messages."""


class PoseGateError(Exception):
    """Base class for expected scientific or input validation failures."""


class ConfigurationError(PoseGateError):
    """The frozen scientific configuration is incomplete or inconsistent."""


class InputValidationError(PoseGateError):
    """Topology, trajectory, selections, timing, or PBC data are invalid."""


class CheckpointNotReady(PoseGateError):
    """The trajectory does not yet contain the configured coordinate prefix."""


class OutcomeEmbargoError(PoseGateError):
    """A genuine shadow record cannot be created after late data are visible."""


class ImmutableRecordError(PoseGateError):
    """An immutable prediction record would be overwritten or is invalid."""
