# src/dna_repeat/error.py
    # Error handling cases

class InvalidFASTAError(Exception):
    exit_code = 1
    def __init__(self, message: str, details: str | None = None) -> None:
        self.message = message
        self.details = details
        super().__init__(message)

class EmptySequenceError(Exception):
    def __init__(self, rec_id: str) -> None:
        full_message = f'{rec_id} : sequence is empty.'
        super().__init__(full_message)

class InvalidSequenceError(Exception):
    def __init__(self, rec_id: str, invalid_chars: list[str]) -> None:
        invalid_chars_str = ', '.join(invalid_chars)
        full_message = f'{rec_id} has invalid bases: {invalid_chars_str}'
        super().__init__(full_message)

class InvalidKmerError(Exception):
    def __init__(self, rec_id: str) -> None:
        full_message = f'{rec_id} : k-mer length is longer than the query sequence.'
        super().__init__(full_message)