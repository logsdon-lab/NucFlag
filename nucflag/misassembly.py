from enum import StrEnum, auto


class Misassembly(StrEnum):
    COLLAPSE_OTHER = auto()
    COLLAPSE_VAR = auto()
    COLLAPSE = auto()
    MISJOIN = auto()
    FALSE_DUPE = auto()

    def as_color(self) -> str:
        match self:
            case self.COLLAPSE_VAR:
                return "blue"
            case self.COLLAPSE:
                return "green"
            case self.MISJOIN:
                return "orange"
            case self.FALSE_DUPE:
                return "purple"
            case self.COLLAPSE_OTHER:
                return "red"
            case _:
                raise ValueError(f"Unreachable {self}")
