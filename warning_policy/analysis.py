"""整理统一基准的水位范围与候选流域关系。"""
from collections import defaultdict


def water_ranges(tables):
    datums = {(row["gauge_id"], row["datum_version"]): row["offset_m"] for row in tables["datums"]}
    values = defaultdict(list)
    unsupported = []
    for row in tables["levels"]:
        key = (row["gauge_id"], row["datum_version"])
        if key not in datums or row["level_m"] is None:
            unsupported.append(row["record_id"])
            continue
        values[row["gauge_id"]].append(row["level_m"] + datums[key])
    return {"gauges": {key: {"minimum_m": min(items), "maximum_m": max(items), "count": len(items)}
                       for key, items in sorted(values.items())}, "unsupported_records": unsupported}


def summarize(tables):
    return {**water_ranges(tables), "candidate_links": len(tables["network"]),
            "scope": "统一基准的水位范围，不评估时滞关联"}
