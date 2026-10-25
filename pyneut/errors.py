"""Errors raised while reading a model."""


class InputError(ValueError):
    """A model input is invalid.

    Carries every problem found, each prefixed with where it is in the input
    (e.g. ``cells.shield.material``), so one run of ``pyneut check`` reports
    all of them rather than the first.
    """

    def __init__(self, problems):
        if isinstance(problems, str):
            problems = [problems]
        self.problems = list(problems)
        super().__init__("\n".join(self.problems))


class ProblemList:
    """Collects input problems and raises them together."""

    def __init__(self):
        self.items = []

    def add(self, where, message):
        self.items.append(f"{where}: {message}" if where else message)

    def raise_if_any(self):
        if self.items:
            raise InputError(self.items)
