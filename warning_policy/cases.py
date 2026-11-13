"""按事件保存预警与影响标签。"""

from dataclasses import dataclass


@dataclass(frozen=True)
class EventCase:
    event_id: str
    score: float
    observed_impact: bool | None


def classify(case: EventCase, threshold: float) -> str:
    if case.observed_impact is None:
        return "unknown"
    return "alert" if case.score >= threshold else "no_alert"
