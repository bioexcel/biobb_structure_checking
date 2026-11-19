"""
 Exceptions raised by biobb_structure_checking.

 All of them derive from StructureCheckingError, so library users (Notebooks,
 biobb wrappers, scripts) can catch a single class. The command line entry
 point (check_structure) converts them into a message on stderr and an exit
 code (StructureCheckingError.exit_code).
"""


class StructureCheckingError(Exception):
    """ Base class for all biobb_structure_checking errors"""
    exit_code = 1

    def __init__(self, message=''):
        super().__init__(message)
        self.message = message

    def __str__(self):
        return str(self.message)


class StructureLoadError(StructureCheckingError):
    """ Raised when the input structure cannot be fetched, downloaded or parsed"""


class AtomLimitError(StructureLoadError):
    """ Raised when the input structure exceeds the allowed number of atoms"""


class UnsupportedFormatError(StructureLoadError, ValueError):
    """ Raised when a structure file format is not supported"""


class SequenceLoadError(StructureCheckingError):
    """ Raised when a FASTA sequence cannot be retrieved or loaded"""


class LibraryLoadError(StructureCheckingError):
    """ Raised when an internal data, residue or forcefield library cannot be loaded"""
    exit_code = 2


class TimeLimitError(StructureCheckingError):
    """ Raised when the allowed execution time is exceeded"""


class CommandNotFoundError(StructureCheckingError):
    """ Raised when an unknown command is requested"""


class CommandListError(StructureCheckingError):
    """ Raised when a command list is not provided or cannot be read"""


class MutationError(StructureCheckingError):
    """ Raised when requested mutations cannot be performed"""


class RenumberError(StructureCheckingError):
    """ Raised on wrong renumbering rules or residue collisions when renumbering"""


class ChainSelectionError(StructureCheckingError):
    """ Raised when a chain selection would leave an empty structure"""


class BuildCoordsError(StructureCheckingError):
    """ Raised when new atom coordinates cannot be built from the residue library"""


class ModellerError(StructureCheckingError):
    """ Raised when Modeller is not available or cannot be set up"""


# Errors that must abort the whole run. StructureChecking._run_method reports
# any other error raised by a command's fix and continues with the next one.
FATAL_ERRORS = (
    StructureLoadError,
    SequenceLoadError,
    LibraryLoadError,
    TimeLimitError,
    CommandNotFoundError,
    CommandListError,
    MutationError,
    RenumberError,
    ChainSelectionError,
    BuildCoordsError,
    ModellerError,
)
