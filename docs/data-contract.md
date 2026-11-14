# 领域资料与供数约定

## 数据含义

rain_mm 为 at 之前 interval_minutes 内累积量，水位 level_m 是瞬时值。水位基准通过 datum_version 及 offset_m 记录，统一基准值为 level_m+offset_m；缺失基准不允许猜测。network 只规定候选流域关系，不证明统计联系。process_id 用于整体留出检验，同过程不同站不能随机拆入训练和验证。

## 现行摘要

`warning_policy/analysis.py` 生成现行资料摘要。源码和运行结果中的 scope 说明其统计边界，不能把未提供的推断当成已有能力。

## 滞后关联研究

`warning_policy/lag_association.py` 提供独立的时滞关联研究，入口：

```sh
python3 -m warning_policy.runner --study lag-association --output /tmp/lag-study.json
```

它与现行摘要分开调用，`--study` 缺省时仍只输出原摘要。研究把同一 `process_id` 内的雨量强度与统一基准水位配成对齐序列，按 `network` 中 `upstream` 的雨量站组成候选站组，在每个水位站上搜索时滞档位；站组与时滞一起搜索，因此用过程内循环平移的族内最大零分布给出搜索校正 p 值。选型与验证按完整过程划分，同一过程不会同时出现在两侧。

解释边界：

- 结果是被校验证过的统计关联证据，不证明雨量导致水位变化，也不代表洪水预测精度。
- 有效时滞范围受雨量采样间隔与时滞档宽限制，是候选区间，不是精确传播时间。
- 雨量窗口按时间重叠加权；覆盖不足的水位时刻直接不配对，不插值成精确观测。
- 不同 `datum_version` 之间不桥接；缺少基准的水位记录计入 `level_excluded`。
- 同一场雨下高度相关的上游雨量站只作为冗余证据；选组采用“相关不显著下降时取最少站数”的简约规则，并单列等价候选。

## 供数接口

`GET /manifest` 返回 generation、tables（每表 count 与 sha256）及 metadata。`GET /tables/<name>?generation=<generation>` 返回同代 rows。供数换代时旧 generation 返回 HTTP 409，不允许将新代记录接入旧分析；不存在的表返回 404。表内容为稳定键唯一的完整同代记录，不是增量事件流。客户端对清单摘要进行核验。摘要使用 UTF-8、排序对象键、紧凑 JSON 的 SHA256，数组顺序保留供数端顺序。

字段与主键声明位于 `warning_policy/schema.json`，合成资料位于 `fixtures/dataset.json`。`metadata.cutoff` 是供数展示截点，不自动为分析选择训练截止点。服务可将整个 dataset 原子替换为新一代，但不在同一 generation 下更改内容。

## 结果保存

`reports.py` 保存现有摘要及来源表摘要，report_id 是内容摘要，不是随机编号。相同内容生成相同编号，临时文件写入完成后原子替换目标文件；历史内容是否保留由调用者选择不同目标路径。服务只监听本机，不接入生产系统。
