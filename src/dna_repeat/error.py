# src/dna_repeat/error.py
    # Error handling cases

class InvalidFASTAError(Exception):
    exit_code = 1
    def __init__(self, message: str, details: str | None = None) -> None:
        self.message = message
        self.details = details
        super().__init__(message)