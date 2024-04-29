from enum import StrEnum


class Misassembly(StrEnum):
    COLLAPSE_VAR = "COLLAPSE_VAR"
    COLLAPSE = "COLLAPSE"
    MISJOIN = "MISJOIN"
    GAP = "GAP"
    FALSE_DUP = "FALSE_DUP"

    def as_color(self) -> str:
        match self:
            case self.COLLAPSE_VAR:
                return "blue"
            case self.COLLAPSE:
                return "green"
            case self.MISJOIN:
                return "orange"
            case self.GAP:
                return "gray"
            case self.FALSE_DUP:
                return "purple"
            case _:
                raise ValueError(f"Invalid color {self}")
