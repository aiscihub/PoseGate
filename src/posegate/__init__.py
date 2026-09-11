"""PoseGate-MD public package."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("posegate-md")
except PackageNotFoundError:
    __version__ = "0.2.0"

__all__ = ["__version__"]
