"""Domain exceptions with actionable command-line messages."""


class TxIDError(Exception):
    """Base class for expected, user-facing TxID failures."""


class AnnotationParseError(TxIDError):
    """An annotation cannot be converted to the common transcript model."""


class ReferenceError(TxIDError):
    """Reference metadata or a contig mapping is invalid."""


class IdentityCollisionError(TxIDError):
    """A truncated public digest maps to more than one full digest."""


class RegistryError(TxIDError):
    """A registry operation failed validation."""

