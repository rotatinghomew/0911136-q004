"""雨量与河道水位的滞后关联证据研究。

本模块只输出统计关联证据：用于判断某个上游候选站组与水位变化之间是否
存在值得进一步采用的滞后线索，不证明因果，也不提供洪水预测精度。

关键边界：
* 雨量记录是 at 之前 interval_minutes 内的累积量，先换算为强度，再用时间
  重叠加权窗口取值；不对雨量或水位做逐点插值，覆盖不足的窗口直接不配对。
* 水位按 datum_version 分段；不同基准之间不桥接，缺少基准的水位记录显式排除。
* 只使用同一 process_id 内的配对；选型与验证按完整过程划分，同一过程不会
  同时出现在两侧。
* 站组与时滞同时在选型过程上搜索，因此用过程内循环平移的蒙特卡洛零分布
  校正反复搜索带来的偶然高分；验证阶段使用固定选型结果单独检验。
* 上游雨量站之间通常高度相关。报告站间相关与边际增益，并优先选择能解释
  数据的更小站组，避免把同一场雨的重复测点当作独立证据。
"""
from __future__ import annotations

import itertools
import math
import random
from collections import defaultdict
from datetime import datetime, timedelta, timezone

DEFAULTS = {
    # 时滞搜索网格：60 分钟一档，0–360 分钟。
    "bin_minutes": 60,
    "max_lag_minutes": 360,
    # 一个水位时刻的雨量窗口至少要有这么大比例被实测记录覆盖才允许配对。
    "min_window_coverage": 0.5,
    # 单个过程-基准段至少这么多个有效水位配对才计算相关。
    "min_pairs_per_segment": 6,
    # 选型、验证各自至少需要的过程数。
    "min_processes": 4,
    # 一个时滞档位上，有效配对覆盖全部水位时刻的最低比例。
    "min_support": 0.5,
    # 申报“存在候选关联”的最低过程内合并相关。
    "min_peak_r": 0.5,
    # 峰值相关至少要比零滞后档位高这么多，才认为滞后本身有解释力。
    "min_lag_gain": 0.10,
    # 有效时滞范围向峰值两侧扩展时允许的相关下降（绝对量）。
    "range_margin": 0.10,
    # 站组选择的简约余量：更大站组必须至少提高这么多相关才值得保留。
    "min_gain": 0.02,
    "max_group_size": 3,
    "max_candidate_groups": 20,
    # 零分布重抽样次数；p 值分辨率约为 1/(permutations+1)。
    "permutations": 300,
    "alpha": 0.05,
    # 站间相关达到该值时提示上游站彼此冗余。
    "redundant_station_r": 0.90,
    # 验证稳定性判据。
    "validation_positive_share": 0.8,
    "validation_r_drop": 0.15,
    "r_bound": 0.999999,
    "constant_epsilon": 1e-12,
    "seed": 20250909,
}

SCOPE = (
    "过程内雨量强度与统一基准水位的滞后关联统计证据；已校正时滞与站组反复搜索，"
    "并在未参与选型的完整过程上复核。该结果不证明雨量导致水位变化，也不代表"
    "洪水预测精度，不能替代河道巡查与本地经验判断。"
)

LIMITATIONS = [
    "相关不等于因果；共同天气过程、流域蓄水与调度都可能同时影响雨量和水位。",
    "有效时滞范围是候选区间，分辨率受雨量记录间隔与时滞档宽限制，不是精确传播时间。",
    "零分布通过过程内循环平移构造，只能校正所选搜索空间内的偶然高分，不能消除未观测混杂。",
    "验证过程数量有限时，复核只能说明“在留出过程上没有明显反转”，不构成总体泛化保证。",
    "搜索校正 p 值的下限由重抽样次数决定；达到下限只说明未观测到更极端的零分布高分。",
    "站点不足、基准缺失或采样覆盖过低的区间会显式排除，不通过插值补充。",
]

NOT_CLAIMED = [
    "不声称雨量关联已经证明因果。",
    "不声称具备洪水预测精度或可替代预警发布。",
    "不把未配对区间插值成精确观测。",
    "不把同一场雨下彼此高度相关的多个雨量站当作独立证据。",
]

METHOD = [
    "雨量累积量按 interval_minutes 换算为强度，并在过程内标准化，消除跨过程的共同季节水平。",
    "在每个水位时刻回看 [lag, lag+档宽) 分钟，用雨量记录的时间重叠加权取窗口强度；覆盖不足则丢弃该时刻。",
    "水位按 datum_version 分段并统一加 offset_m；不同基准之间不拼接、不桥接。",
    "只在同一 process_id 内把窗口雨量与水位配对，再按基准段计算相关，用 Fisher-z 合并。",
    "在选型过程上同时搜索候选站组与时滞；用过程内循环平移的族内最大零分布校正反复搜索。",
    "按站数最少且相关不显著下降的简约规则给出候选站组，并报告站间相关与冗余提示。",
    "把固定站组与固定时滞档拿到未参与选型的完整过程上复核，检查方向和强度是否稳定。",
]


# ---------------------------------------------------------------- 基础工具


def _label(value):
    """统一把标识转为稳定的字符串，便于排序和输出。"""
    return str(value)


def _text(value):
    return None if value is None else str(value)


def _parse_time(value):
    if isinstance(value, datetime):
        moment = value
    elif isinstance(value, str):
        moment = datetime.fromisoformat(value)
    else:
        raise ValueError("时间字段不是 ISO 时间")
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc)


def _number(value, positive=False):
    if value is None or isinstance(value, bool):
        raise ValueError("不是数值")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("不是有限数值")
    if positive and result <= 0:
        raise ValueError("不是正数")
    return result


def _round(value, digits=4):
    if value is None:
        return None
    return round(float(value), digits)


def _iso(moment):
    return None if moment is None else moment.isoformat()


def _pearson(xs, ys, minimum, epsilon):
    n = len(xs)
    if n < minimum:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx <= epsilon or syy <= epsilon:
        return None
    value = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)
    return max(-1.0, min(1.0, value))


def _fisher(items, minimum, bound):
    """按 (r, n) 的 Fisher-z 权重合并；只接受样本数达标且非常数的段。"""
    kept = []
    for value, size in items:
        if value is None or size < minimum or not (-1.0 < value < 1.0):
            continue
        kept.append((value, size))
    weight = sum(size - 3 for _, size in kept)
    if weight <= 0:
        return None, 0, []
    z_value = sum((size - 3) * math.atanh(max(-bound, min(bound, value)))
                  for value, size in kept) / weight
    return math.tanh(z_value), weight, kept


def _median(values):
    ordered = sorted(values)
    if not ordered:
        return None
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


# ---------------------------------------------------------------- 数据准备


def _datum_summary(segments):
    """按 (水位站, 基准版本) 汇总分段，避免把每一场过程都铺开。"""
    grouped = defaultdict(list)
    for item in segments:
        grouped[(item["gauge_id"], item["datum_version"])].append(item)
    summary = []
    for (gauge, version), items in sorted(grouped.items()):
        items.sort(key=lambda item: item["start"])
        summary.append({
            "gauge_id": gauge,
            "datum_version": version,
            "segments": len(items),
            "processes": sorted({item["process_id"] for item in items}),
            "records": sum(item["records"] for item in items),
            "start": items[0]["start"],
            "end": max(item["end"] for item in items),
        })
    return summary


def _empty_exclusions():
    return defaultdict(int)


def _prepare(tables, config):
    """把两张观测表整理成可按过程配对的结构，并记录所有排除原因。"""
    rain_excluded = _empty_exclusions()
    level_excluded = _empty_exclusions()
    rain_total = len(tables.get("rain", []))
    level_total = len(tables.get("levels", []))
    rain_records = defaultdict(list)
    level_records = defaultdict(list)
    station_intervals = defaultdict(list)
    station_used = defaultdict(int)
    datum_offsets = {}

    for row in tables.get("datums", []):
        gauge = _label(row.get("gauge_id"))
        version = _label(row.get("datum_version"))
        try:
            offset = _number(row.get("offset_m"))
        except (TypeError, ValueError):
            continue
        datum_offsets[(gauge, version)] = offset

    for row in tables.get("rain", []):
        process = row.get("process_id")
        station = row.get("station_id")
        if process is None:
            rain_excluded["missing_process"] += 1
            continue
        if station is None:
            rain_excluded["missing_station"] += 1
            continue
        try:
            end = _parse_time(row.get("at"))
        except (TypeError, ValueError):
            rain_excluded["invalid_time"] += 1
            continue
        try:
            amount = _number(row.get("rain_mm"))
        except (TypeError, ValueError):
            rain_excluded["missing_or_invalid_amount"] += 1
            continue
        try:
            interval = _number(row.get("interval_minutes"), positive=True)
        except (TypeError, ValueError):
            rain_excluded["missing_or_invalid_interval"] += 1
            continue
        record = {
            "at": end,
            "end": end,
            "start": end - timedelta(minutes=interval),
            "interval": interval,
            "amount": amount,
            "rate": amount * 60.0 / interval,
            "record_id": _text(row.get("record_id")),
        }
        rain_records[(_label(process), _label(station))].append(record)
        station_intervals[_label(station)].append(interval)

    for row in tables.get("levels", []):
        process = row.get("process_id")
        gauge = row.get("gauge_id")
        if process is None:
            level_excluded["missing_process"] += 1
            continue
        if gauge is None:
            level_excluded["missing_gauge"] += 1
            continue
        try:
            moment = _parse_time(row.get("at"))
        except (TypeError, ValueError):
            level_excluded["invalid_time"] += 1
            continue
        try:
            level = _number(row.get("level_m"))
        except (TypeError, ValueError):
            level_excluded["missing_or_invalid_level"] += 1
            continue
        version = _label(row.get("datum_version"))
        offset = datum_offsets.get((_label(gauge), version))
        if offset is None:
            level_excluded["missing_datum_offset"] += 1
            continue
        level_records[(_label(process), _label(gauge))].append({
            "at": moment,
            "value": level + offset,
            "datum": version,
            "record_id": _text(row.get("record_id")),
        })

    rain_series = {}
    rain_notes = defaultdict(list)
    for key, records in rain_records.items():
        records.sort(key=lambda item: (item["start"], item["at"], item["record_id"] or ""))
        moments = [item["at"] for item in records]
        if len(set(moments)) != len(moments):
            rain_notes[key].append("ambiguous_timestamps")
            rain_excluded["ambiguous_timestamps"] += len(records)
            continue
        overlap = any(records[index]["start"] < records[index - 1]["end"]
                      for index in range(1, len(records)))
        if overlap:
            # 区间累积量互相重叠时无法判断哪一段才是实测值，整条序列不参与配对。
            rain_notes[key].append("overlapping_intervals")
            rain_excluded["overlapping_intervals"] += len(records)
            continue
        weight = sum(item["interval"] for item in records)
        mean = sum(item["rate"] * item["interval"] for item in records) / weight
        variance = sum(item["interval"] * (item["rate"] - mean) ** 2 for item in records) / weight
        spread = math.sqrt(max(variance, 0.0))
        if spread <= config["constant_epsilon"]:
            rain_notes[key].append("constant_intensity")
            rain_excluded["constant_intensity"] += len(records)
            continue
        for item in records:
            item["z"] = (item["rate"] - mean) / spread
        rain_series[key] = records
        station_used[key[1]] += len(records)

    level_series = {}
    level_notes = defaultdict(list)
    datum_segments = []
    for key, records in level_records.items():
        records.sort(key=lambda item: (item["at"], item["record_id"] or ""))
        grouped = defaultdict(list)
        for item in records:
            grouped[item["at"]].append(item)
        unique = []
        ambiguous = 0
        for moment in sorted(grouped):
            values = grouped[moment]
            if len(values) > 1:
                ambiguous += len(values)
                continue
            unique.append(values[0])
        if ambiguous:
            level_notes[key].append("ambiguous_timestamps")
            level_excluded["ambiguous_timestamps"] += ambiguous
        if not unique:
            continue
        segments = []
        current = [unique[0]]
        for item in unique[1:]:
            if item["datum"] != current[-1]["datum"]:
                segments.append(current)
                current = [item]
            else:
                current.append(item)
        segments.append(current)
        prepared = []
        for segment in segments:
            datum_segments.append({
                "process_id": key[0],
                "gauge_id": key[1],
                "datum_version": segment[0]["datum"],
                "start": _iso(segment[0]["at"]),
                "end": _iso(segment[-1]["at"]),
                "records": len(segment),
            })
            prepared.append({
                "datum": segment[0]["datum"],
                "points": [{"at": item["at"], "value": item["value"]} for item in segment],
            })
        level_series[key] = prepared

    return {
        "rain_series": rain_series,
        "rain_notes": {key: sorted(values) for key, values in sorted(rain_notes.items())},
        "level_series": level_series,
        "level_notes": {key: sorted(values) for key, values in sorted(level_notes.items())},
        "rain_excluded": dict(sorted(rain_excluded.items())),
        "level_excluded": dict(sorted(level_excluded.items())),
        "rain_total": rain_total,
        "level_total": level_total,
        "station_intervals": {
            key: {
                "records": len(values),
                "min_minutes": min(values),
                "median_minutes": _median(values),
                "max_minutes": max(values),
                "distinct_minutes": sorted(set(values)),
            } for key, values in sorted(station_intervals.items())
        },
        "station_used": dict(sorted(station_used.items())),
        "datum_segments": _datum_summary(datum_segments),
    }


# ---------------------------------------------------------------- 配对与前向取值


def _window_value(records, moment, low, high, config):
    """取 [moment-high, moment-low] 内按时间重叠加权的雨量强度。

    记录本身就是区间累积量，因此按完全落在窗口内的实测时段加权；
    窗口被实测覆盖的比例不足时返回 None，绝不插值成“精确配对”。
    """
    window_start = moment - timedelta(minutes=high)
    window_end = moment - timedelta(minutes=low)
    numerator = 0.0
    covered = 0.0
    for record in records:
        start = record["start"] if record["start"] > window_start else window_start
        end = record["end"] if record["end"] < window_end else window_end
        if end <= start:
            continue
        overlap = (end - start).total_seconds() / 60.0
        numerator += record["z"] * overlap
        covered += overlap
    if covered <= 0:
        return None, 0.0
    span = max(high - low, 1e-9)
    fraction = covered / span
    if fraction < config["min_window_coverage"]:
        return None, fraction
    return numerator / covered, fraction


def _bin_edges(config):
    """最远只看 max_lag_minutes，最后一个档位完整落在该上限内。"""
    step = int(config["bin_minutes"])
    limit = int(config["max_lag_minutes"])
    return [(low, low + step) for low in range(0, limit, step)]


def _align_segments(prepared, gauge, group, low, high, processes, config):
    """按 (过程, 基准段) 构造对齐序列；覆盖不足的水位时刻整体不配对。"""
    aligned = []
    for process in processes:
        segments = prepared["level_series"].get((process, gauge), [])
        for index, segment in enumerate(segments):
            xs = []
            ys = []
            indices = []
            total = len(segment["points"])
            for position, point in enumerate(segment["points"]):
                values = []
                for station in group:
                    records = prepared["rain_series"].get((process, station))
                    if not records:
                        values = []
                        break
                    value, _ = _window_value(records, point["at"], low, high, config)
                    if value is None:
                        values = []
                        break
                    values.append(value)
                if values:
                    xs.append(sum(values) / len(values))
                    ys.append(point["value"])
                    indices.append(position)
            aligned.append({
                "process_id": process,
                "datum_version": segment["datum"],
                "segment_index": index,
                "xs": xs,
                "ys": ys,
                "indices": indices,
                "full_ys": [point["value"] for point in segment["points"]],
                "total": total,
                "pairs": len(xs),
            })
    return aligned


def _segment_correlation(item, config, ys=None):
    values = item["ys"] if ys is None else ys
    return _pearson(item["xs"], values, config["min_pairs_per_segment"],
                    config["constant_epsilon"])


def _aggregate(aligned, config):
    """把各段相关按 Fisher-z 合并，并给出有效配对覆盖率。"""
    items = []
    supported = 0
    total = 0
    used_processes = set()
    for item in aligned:
        total += item["total"]
        if item["pairs"] < config["min_pairs_per_segment"]:
            continue
        ys = item["ys"]
        correlation = _segment_correlation(item, config, ys)
        if correlation is None:
            continue
        items.append((correlation, item["pairs"]))
        supported += item["pairs"]
        used_processes.add(item["process_id"])
    combined, weight, kept = _fisher(items, config["min_pairs_per_segment"], config["r_bound"])
    support = supported / total if total else 0.0
    return combined, weight, len(kept), support, sorted(used_processes)


def _profiles(prepared, gauge, groups, processes, edges, config):
    """为每个候选站组在每个时滞档位上计算选型期相关与覆盖。"""
    result = {}
    for group in groups:
        for low, high in edges:
            aligned = _align_segments(prepared, gauge, group, low, high, processes, config)
            combined, weight, segments, support, used = _aggregate(aligned, config)
            result[(group, low, high)] = {
                "aligned": aligned,
                "r": combined,
                "weight": weight,
                "segments": segments,
                "support": support,
                "processes": used,
            }
    return result


def _admissible(entry, config):
    return (entry["r"] is not None
            and entry["support"] >= config["min_support"]
            and len(entry["processes"]) >= config["min_processes"])


def _peak(profiles, group, edges, config):
    best = None
    for low, high in edges:
        entry = profiles[(group, low, high)]
        if not _admissible(entry, config):
            continue
        candidate = (entry["r"], -low, low, high, entry)
        if best is None or candidate[:4] > best[:4]:
            best = candidate
    return best


def _effective_range(profiles, group, edges, peak_low, peak_r, config):
    qualifying = []
    for low, high in edges:
        entry = profiles[(group, low, high)]
        if _admissible(entry, config) and entry["r"] >= peak_r - config["range_margin"]:
            qualifying.append((low, high))
    if not qualifying:
        return None
    qualifying.sort()
    start = next(index for index, item in enumerate(qualifying) if item[0] == peak_low)
    left = start
    while left > 0 and qualifying[left - 1][1] == qualifying[left][0]:
        left -= 1
    right = start
    while right + 1 < len(qualifying) and qualifying[right + 1][0] == qualifying[right][1]:
        right += 1
    return {
        "low_minutes": qualifying[left][0],
        "high_minutes": qualifying[right][1],
        "bins": [{"low_minutes": low, "high_minutes": high,
                  "r": _round(profiles[(group, low, high)]["r"])} for low, high in qualifying[left:right + 1]],
    }


# ---------------------------------------------------------------- 站间相关


def _station_pair_correlation(prepared, gauge, stations, processes, config):
    """在选型过程上测算两站雨量强度的相关（诊断用，不作为独立样本）。

    各站采样时刻不同，因此不按时间戳硬配对；改用两站中较长的采样间隔作为
    共同诊断窗口，在水位观测时刻上分别取重叠加权值，窗口覆盖不足的时刻
    直接丢弃。逐过程计算相关后取中位数。
    """
    left, right = stations
    spans = []
    for station in stations:
        intervals = [record["interval"] for process in processes
                     for record in prepared["rain_series"].get((process, station), [])]
        if not intervals:
            return None, 0
        spans.append(_median(intervals))
    window = max(spans)
    per_process = []
    for process in processes:
        left_records = prepared["rain_series"].get((process, left))
        right_records = prepared["rain_series"].get((process, right))
        if not left_records or not right_records:
            continue
        xs = []
        ys = []
        for segment in prepared["level_series"].get((process, gauge), []):
            for point in segment["points"]:
                value_left, _ = _window_value(left_records, point["at"], 0.0, window, config)
                value_right, _ = _window_value(right_records, point["at"], 0.0, window, config)
                if value_left is not None and value_right is not None:
                    xs.append(value_left)
                    ys.append(value_right)
        value = _pearson(xs, ys, config["min_pairs_per_segment"], config["constant_epsilon"])
        if value is not None:
            per_process.append(value)
    if not per_process:
        return None, 0
    return _median(per_process), len(per_process)


def _station_correlation_matrix(prepared, gauge, stations, processes, config):
    matrix = {}
    for left, right in itertools.combinations(sorted(stations), 2):
        value, count = _station_pair_correlation(prepared, gauge, (left, right), processes, config)
        matrix[(left, right)] = {"r": _round(value), "processes": count}
    return matrix


# ---------------------------------------------------------------- 零分布


def _rotation_rng(seed, gauge):
    material = f"{seed}|{gauge}"
    return random.Random(material)


def _null_family_max(label, profiles, groups, edges, processes, config, permutations):
    """过程内循环平移水位序列，保留过程内形状，打散与雨量的时间对应。

    每个重抽样都重新搜索全部站组与时滞，取族内最大相关，才能反映
    “反复搜索时滞和站组”这一实际选择过程带来的偶然高分。对齐关系在
    选型期已经固定，重抽样只旋转水位值，不重新选择配对，避免改变覆盖结构。
    """
    rng = _rotation_rng(config["seed"], label)
    first_low, first_high = edges[0]
    template = profiles[(groups[0], first_low, first_high)]["aligned"]
    sizes = {(item["process_id"], item["segment_index"]): len(item["full_ys"])
             for item in template}
    maxima = []
    for _ in range(permutations):
        # 同一 (过程, 基准段) 只用一个平移量，保持组间与档位间共用同一条零水位序列。
        rotations = {key: (rng.randrange(1, size) if size > 1 else 0)
                     for key, size in sizes.items()}
        limit = -2.0
        for group in groups:
            for low, high in edges:
                shifted = []
                for item in profiles[(group, low, high)]["aligned"]:
                    shift = rotations[(item["process_id"], item["segment_index"])]
                    full = item["full_ys"]
                    moved = [full[(index + shift) % len(full)] for index in item["indices"]]
                    shifted.append({**item, "ys": moved})
                combined, _, _, support, used = _aggregate(shifted, config)
                if combined is None or support < config["min_support"]:
                    continue
                if len(used) < config["min_processes"]:
                    continue
                if combined > limit:
                    limit = combined
        maxima.append(limit)
    maxima.sort()
    return maxima


def _p_value(observed, maxima):
    if observed is None or not maxima:
        return None
    exceed = sum(1 for value in maxima if value >= observed)
    return (exceed + 1) / (len(maxima) + 1)


# ---------------------------------------------------------------- 选型


def _split_processes(prepared, gauges):
    """按过程起始时间排序后整段留出；绝不把同一过程的站点拆到两侧。"""
    starts = {}
    for (process, _station), records in prepared["rain_series"].items():
        if records:
            starts[process] = min(starts.get(process, records[0]["at"]), records[0]["at"])
    for (process, _gauge), segments in prepared["level_series"].items():
        if segments:
            first = segments[0]["points"][0]["at"]
            starts[process] = min(starts.get(process, first), first)
    ordered = sorted(starts, key=lambda name: (starts[name], name))
    usable = [name for name in ordered
              if any((name, gauge) in prepared["level_series"] for gauge in gauges)]
    return usable


def _candidate_groups(stations, config):
    groups = []
    for size in range(1, min(config["max_group_size"], len(stations)) + 1):
        for combination in itertools.combinations(sorted(stations), size):
            groups.append(combination)
    groups.sort(key=lambda group: (len(group), group))
    return groups[:config["max_candidate_groups"]], len(groups) > config["max_candidate_groups"]


def _rank_groups(profiles, groups, edges, catchments, config):
    ranked = []
    for group in groups:
        found = _peak(profiles, group, edges, config)
        if found is None:
            continue
        r_value, _neg_low, low, high, entry = found
        ranked.append({
            "group": group,
            "r": r_value,
            "low_minutes": low,
            "high_minutes": high,
            "support": entry["support"],
            "processes": entry["processes"],
            "catchments": sorted({catchments[station] for station in group
                                  if catchments.get(station) is not None}),
        })
    ranked.sort(key=lambda item: (-item["r"], len(item["group"]), item["group"], item["low_minutes"]))
    return ranked


def _select(ranked, config):
    """选出一个可解释的站组：更大站组必须带来超过余量的边际增益。

    先取全部候选中的最高相关，再在“相关不超过最高值减去 min_gain”的范围内
    选择站数最少者；站数相同时优先覆盖更多水位时刻、跨更多支流、相关更高
    的候选。这样不会因为多加一个高度相关的雨量站而把同一场雨的重复测点
    当成独立证据。
    """
    if not ranked:
        return None
    best_r = max(item["r"] for item in ranked)
    eligible = [item for item in ranked if item["r"] >= best_r - config["min_gain"]]
    eligible.sort(key=lambda item: (len(item["group"]), -item["support"],
                                    -len(item["catchments"]), -item["r"],
                                    item["group"], item["low_minutes"]))
    selected = eligible[0]
    alternatives = [item for item in eligible
                    if len(item["group"]) == len(selected["group"])
                    and item["group"] != selected["group"]]
    return selected, alternatives


def _station_contributions(prepared, gauge, group, low, high, processes, config):
    """单站分别的合并相关，用于说明组内每站的独立贡献（不宣称独立性）。"""
    result = {}
    for station in group:
        aligned = _align_segments(prepared, gauge, (station,), low, high, processes, config)
        combined, _, _, support, used = _aggregate(aligned, config)
        result[station] = {"r": _round(combined), "support": _round(support),
                           "processes": len(used)}
    return result


# ---------------------------------------------------------------- 验证


def _validate(prepared, gauge, selection, processes, config):
    """在未参与选型的完整过程上复核固定站组与固定时滞档。"""
    low = selection["low_minutes"]
    high = selection["high_minutes"]
    aligned = _align_segments(prepared, gauge, selection["group"], low, high, processes, config)
    combined, weight, segments, support, used = _aggregate(aligned, config)
    per_process = []
    for item in aligned:
        if item["pairs"] < config["min_pairs_per_segment"]:
            per_process.append({"process_id": item["process_id"],
                                "datum_version": item["datum_version"],
                                "r": None, "pairs": item["pairs"],
                                "note": "insufficient_pairs"})
            continue
        value = _segment_correlation(item, config)
        per_process.append({"process_id": item["process_id"],
                            "datum_version": item["datum_version"],
                            "r": _round(value), "pairs": item["pairs"],
                            "note": None if value is not None else "constant_level"})
    usable = [item["r"] for item in per_process if item["r"] is not None]
    positive = sum(1 for value in usable if value > 0)
    share = positive / len(usable) if usable else None
    stable = (combined is not None
              and len(usable) >= config["min_processes"]
              and support >= config["min_support"]
              and share is not None and share >= config["validation_positive_share"]
              and combined >= selection["r"] - config["validation_r_drop"])
    return {
        "processes": processes,
        "r": _round(combined),
        "support": _round(support),
        "segments": segments,
        "positive_processes": positive,
        "usable_processes": len(usable),
        "positive_share": _round(share),
        "stable": bool(stable),
        "per_process": per_process,
    }


# ---------------------------------------------------------------- 汇总


def _status(selection, p_value, validation, zero_r, config):
    if selection is None:
        return "no_supported_candidate"
    if selection["r"] < config["min_peak_r"]:
        return "below_peak_threshold"
    if p_value is None or p_value > config["alpha"]:
        return "not_significant_after_search"
    if not validation["stable"]:
        return "unstable_on_held_out_processes"
    if zero_r is not None and selection["r"] - zero_r < config["min_lag_gain"]:
        # 峰值与“水位前最近一档”差不多时，所谓滞后并不比零滞后更有解释力。
        return "not_distinct_from_zero_lag"
    return "supported_candidate"


def _recommendation(status, low, high):
    if status == "supported_candidate":
        return (f"可将上游站组在约 {low}–{high} 分钟前的雨量作为观察线索，"
                "并结合现场巡查复核；该区间是统计候选范围，不是精确传播时间。")
    if status == "unstable_on_held_out_processes":
        return "选型期关联在留出过程上未稳定重现，建议暂不据此调整关注提前量。"
    if status == "not_distinct_from_zero_lag":
        return "峰值与零滞后档位接近，无法据此说明需要提前关注，暂不调整提前量。"
    if status == "not_significant_after_search":
        return "在反复搜索站组与时滞后未超过偶然高分范围，暂不作为时滞线索。"
    if status == "below_peak_threshold":
        return "关联强度不足以形成可用的时滞线索。"
    if status == "insufficient_processes":
        return "可用过程数不足，无法按过程划分选型与验证，暂不给出时滞线索。"
    return "没有同时满足覆盖与过程数要求的候选站组，无法给出时滞范围。"


def _empty_result(gauge, status, selection_processes, validation_processes, note, config):
    """所有拒绝状态共享同一输出骨架，调用方无需判断键是否存在。"""
    return {
        "gauge_id": gauge,
        "status": status,
        "selection_processes": selection_processes,
        "validation_processes": validation_processes,
        "note": note,
        "lag_profile": [],
        "candidate_groups": [],
        "selection": None,
        "top_raw_group": None,
        "candidate_search_truncated": False,
        "zero_lag_reference": None,
        "effective_lag": None,
        "station_correlations": {},
        "redundant_stations": False,
        "search_adjusted_p": None,
        "p_resolution": _round(1.0 / (config["permutations"] + 1)),
        "null_max_median": None,
        "null_max_quantile": None,
        "validation": None,
        "recommendation": _recommendation(status, None, None),
    }


def _gauge_analysis(prepared, gauge, stations, catchments, selection_processes,
                    validation_processes, config):
    edges = _bin_edges(config)
    groups, truncated = _candidate_groups(stations, config)
    if len(selection_processes) < config["min_processes"] or len(validation_processes) < config["min_processes"]:
        return _empty_result(gauge, "insufficient_processes", selection_processes,
                             validation_processes,
                             "选型或验证过程数不足，不做搜索，避免小样本高分。",
                             config), None

    profiles = _profiles(prepared, gauge, groups, selection_processes, edges, config)
    ranked = _rank_groups(profiles, groups, edges, catchments, config)
    top_group = ranked[0]["group"] if ranked else None
    choice = _select(ranked, config)
    selected, alternatives = choice if choice else (None, [])
    if selected is None:
        return _empty_result(gauge, "no_supported_candidate", selection_processes,
                             validation_processes,
                             "没有任何站组在最低覆盖与过程数要求下形成可用相关。",
                             config), None

    maxima = _null_family_max(gauge, profiles, groups, edges, selection_processes, config,
                              config["permutations"])
    p_value = _p_value(selected["r"], maxima)
    validation = _validate(prepared, gauge, selected, validation_processes, config)
    zero_entry = profiles[(selected["group"], edges[0][0], edges[0][1])]
    status = _status(selected, p_value, validation, zero_entry["r"], config)

    effective = _effective_range(profiles, selected["group"], edges,
                                 selected["low_minutes"], selected["r"], config)
    contributions = _station_contributions(prepared, gauge, selected["group"],
                                           selected["low_minutes"], selected["high_minutes"],
                                           selection_processes, config)
    matrix = _station_correlation_matrix(prepared, gauge, stations,
                                         selection_processes, config)
    redundant = any(item["r"] is not None and item["r"] >= config["redundant_station_r"]
                    for item in matrix.values())
    profile = []
    for low, high in edges:
        entry = profiles[(selected["group"], low, high)]
        inside = (effective is not None
                  and effective["low_minutes"] <= low and high <= effective["high_minutes"])
        profile.append({
            "low_minutes": low,
            "high_minutes": high,
            "r": _round(entry["r"]),
            "support": _round(entry["support"]),
            "unpaired_share": _round(1.0 - entry["support"]),
            "segments": entry["segments"],
            "processes": len(entry["processes"]),
            "within_range": bool(inside),
        })
    candidates = []
    for item in ranked[:config["max_candidate_groups"]]:
        candidates.append({
            "stations": list(item["group"]),
            "catchments": item["catchments"],
            "cross_catchment": len(item["catchments"]) > 1,
            "r": _round(item["r"]),
            "lag_minutes": [item["low_minutes"], item["high_minutes"]],
            "support": _round(item["support"]),
            "processes": len(item["processes"]),
            "selected": item["group"] == selected["group"],
        })
    zero_profile = {
        "r": _round(zero_entry["r"]),
        "support": _round(zero_entry["support"]),
        "processes": len(zero_entry["processes"]),
        "note": "0 档位表示水位时刻前的最近一个档宽内即已下雨，滞后证据需明显高于此值",
    }
    low = effective["low_minutes"] if effective else None
    high = effective["high_minutes"] if effective else None
    return {
        "gauge_id": gauge,
        "zero_lag_reference": zero_profile,
        "status": status,
        "selection_processes": selection_processes,
        "validation_processes": validation_processes,
        "selection": {
            "stations": list(selected["group"]),
            "r": _round(selected["r"]),
            "lag_minutes": [selected["low_minutes"], selected["high_minutes"]],
            "support": _round(selected["support"]),
            "processes": len(selected["processes"]),
            "catchments": selected["catchments"],
            "cross_catchment": len(selected["catchments"]) > 1,
            "station_contributions": contributions,
            "equivalent_alternatives": [
                {"stations": list(item["group"]), "r": _round(item["r"]),
                 "lag_minutes": [item["low_minutes"], item["high_minutes"]]}
                for item in alternatives
            ],
        },
        "top_raw_group": list(top_group) if top_group else None,
        "candidate_search_truncated": truncated,
        "effective_lag": None if effective is None else {
            "low_minutes": effective["low_minutes"],
            "high_minutes": effective["high_minutes"],
            "bins": effective["bins"],
            "resolution_minutes": config["bin_minutes"],
            "criterion": (f"以 {config['bin_minutes']} 分钟为档宽，峰值相关下探不超过 "
                          f"{config['range_margin']} 且覆盖达标的连续档位；档宽即分辨率上限"),
        },
        "lag_profile": profile,
        "candidate_groups": candidates,
        "station_correlations": {
            f"{left}|{right}": item for (left, right), item in sorted(matrix.items())
        },
        "redundant_stations": redundant,
        "search_adjusted_p": _round(p_value),
        "p_resolution": _round(1.0 / (config["permutations"] + 1)),
        "null_max_median": _round(_median(maxima)),
        "null_max_quantile": _round(maxima[int(0.95 * (len(maxima) - 1))]) if maxima else None,
        "validation": validation,
        "recommendation": _recommendation(status, low, high),
        "note": ("站组按“相关不超过最高值减去 %.2f 时选站数最少者”的简约规则选择；"
                 "站数相同时优先跨不同支流的组合。站间高度相关时该组只代表冗余证据；"
                 "等价候选另列，不表示可独立复现。"
                 % config["min_gain"]),
    }, maxima


def analyze_lag_association(tables, config=None, selection_share=0.5):
    """生成时滞关联候选站组与有效时滞范围，并在留出过程上复核。

    selection_share 是按过程起始时间排序后用于选型的比例；验证始终使用
    其余完整过程，且同一 process_id 不会跨两侧。
    """
    options = dict(DEFAULTS)
    if config:
        unknown = sorted(set(config) - set(DEFAULTS))
        if unknown:
            raise ValueError("未知的滞后关联参数：" + ", ".join(unknown))
        options.update(config)
    if not 0 < selection_share < 1:
        raise ValueError("选型比例必须在 0 与 1 之间")
    for name in ("bin_minutes", "max_lag_minutes", "min_processes", "permutations",
                 "max_group_size", "max_candidate_groups", "min_pairs_per_segment"):
        if not isinstance(options[name], int) or options[name] < 1:
            raise ValueError("参数必须是正整数：" + name)
    prepared = _prepare(tables, options)
    gauges = sorted({gauge for (_process, gauge) in prepared["level_series"]})
    stations = sorted({station for (_process, station) in prepared["rain_series"]})
    processes = _split_processes(prepared, gauges)
    if len(processes) < 2 * options["min_processes"]:
        # 过程数不足时不做选型，明确输出“无法划分”而不是用少量过程硬搜。
        selection_processes = []
        validation_processes = []
    else:
        split = int(round(len(processes) * selection_share))
        split = max(options["min_processes"], split)
        split = min(split, len(processes) - options["min_processes"])
        selection_processes = processes[:split]
        validation_processes = processes[split:]

    network = defaultdict(list)
    for row in tables.get("network", []):
        gauge = _label(row.get("gauge_id"))
        station = _label(row.get("rain_station"))
        if gauge in gauges and station in stations:
            network[gauge].append({
                "rain_station": station,
                "catchment": _text(row.get("catchment")),
                "upstream": bool(row.get("upstream")),
            })
    results = []
    maxima_by_gauge = {}
    for gauge in gauges:
        links = sorted(network.get(gauge, []), key=lambda item: (item["rain_station"],))
        upstream = sorted({item["rain_station"] for item in links if item["upstream"]})
        candidates = upstream or sorted({item["rain_station"] for item in links})
        if not candidates:
            candidates = stations
        catchments = {item["rain_station"]: item["catchment"] for item in links}
        result, maxima = _gauge_analysis(prepared, gauge, candidates, catchments,
                                         selection_processes, validation_processes, options)
        result["candidate_stations"] = candidates
        result["candidate_basis"] = ("network.upstream" if upstream else
                                     "network.all" if links else "all_stations")
        result["network_links"] = links
        results.append(result)
        if maxima is not None:
            maxima_by_gauge[gauge] = maxima

    return {
        "scope": SCOPE,
        "not_claimed": NOT_CLAIMED,
        "limitations": LIMITATIONS,
        "selection_share": selection_share,
        "test_family": "每个水位站各自一个搜索族；跨水位站的 p 值未做多重比较校正",
        "method": METHOD,
        "config": {key: options[key] for key in sorted(options)},
        "data": {
            "rain_records": prepared["rain_total"],
            "level_records": prepared["level_total"],
            "rain_excluded": prepared["rain_excluded"],
            "level_excluded": prepared["level_excluded"],
            "rain_notes": prepared["rain_notes"],
            "level_notes": prepared["level_notes"],
            "station_intervals": prepared["station_intervals"],
            "station_used": prepared["station_used"],
            "datum_segments": prepared["datum_segments"],
            "processes": processes,
            "selection_processes": selection_processes,
            "validation_processes": validation_processes,
        },
        "gauges": results,
    }
