from enum import StrEnum


class Misassembly(StrEnum):
    COLLAPSE_OTHER = "COLLAPSE_OTHER"
    COLLAPSE_VAR = "COLLAPSE_VAR"
    COLLAPSE = "COLLAPSE"
    MISJOIN = "MISJOIN"
    FALSE_DUP = "FALSE_DUP"
    HET = "HET"

    def as_color(self) -> str:
        match self:
            case self.COLLAPSE_VAR:
                return "blue"
            case self.COLLAPSE:
                return "green"
            case self.MISJOIN:
                return "orange"
            case self.FALSE_DUP:
                return "purple"
            case self.HET:
                return "teal"
            case self.COLLAPSE_OTHER:
                return "red"
            case _:
                raise ValueError(f"Invalid color {self}")
