#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd

# =========================
# 1. 文件路径设置
# =========================
input_file = "merged_result.xlsx"
output_file1 = "merge_with_scores.xlsx"
output_file2 = "merge_updated.xlsx"

# support count 的阈值
# percentile >= 75 表示该模型支持“高破坏性”
SUPPORT_THRESHOLD = 75

# =========================
# 2. 读取Excel
# =========================
if not os.path.exists(input_file):
    raise FileNotFoundError(f"未找到文件: {input_file}")

df = pd.read_excel(input_file)

# 模型列配置
id_col = "Gene_name"
model_cols = ["evo2", "plantCAD2_maize", "esm2", "esm_1V", "saport", "ProtSSN"]

# 检查列是否存在
missing_cols = [col for col in [id_col] + model_cols if col not in df.columns]
if missing_cols:
    raise ValueError(f"Excel 中缺少以下列: {missing_cols}")

# 转成数值，无法转换的设为 NaN
for col in model_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# =========================
# 3. 定义函数：把“越低越有害”的原始分数
#    转成“越高越有害”的 percentile
# =========================
def damaging_percentile(series):
    """
    输入：原始分数，且分数越低越有害
    输出：0-100 的 percentile，越高越有害
    """
    s = series.copy()
    valid = s.dropna()

    # 如果全是缺失
    if len(valid) == 0:
        return pd.Series(np.nan, index=series.index)

    # 如果只有1个有效值
    if len(valid) == 1:
        out = pd.Series(np.nan, index=series.index)
        out.loc[valid.index] = 100.0
        return out

    # rank: 最小值 rank=1（最有害）
    ranks = valid.rank(method="average", ascending=True)

    # 转成百分位：最有害=100，最不有害=0
    percentiles = (len(valid) - ranks) / (len(valid) - 1) * 100

    out = pd.Series(np.nan, index=series.index)
    out.loc[valid.index] = percentiles
    return out

# =========================
# 4. 计算各模型 percentile
# =========================
percentile_cols = []
for col in model_cols:
    pct_col = f"{col}_percentile"
    df[pct_col] = damaging_percentile(df[col])
    percentile_cols.append(pct_col)

# =========================
# 5. 计算 consensus score
#    取6个模型 percentile 的平均值
# =========================
df["consensus score"] = df[percentile_cols].mean(axis=1, skipna=True)

# =========================
# 6. 计算 support count
#    percentile >= SUPPORT_THRESHOLD 记为支持
# =========================
df["support count"] = (df[percentile_cols] >= SUPPORT_THRESHOLD).sum(axis=1)

# =========================
# 7. 计算 consistency SD
#    百分位标准差，越小表示模型越一致
# =========================
df["consistency SD"] = df[percentile_cols].std(axis=1, skipna=True, ddof=0)

# 如果只有1个模型有值，std会是0；如果全缺失则保留NaN
available_counts = df[percentile_cols].notna().sum(axis=1)
df.loc[available_counts <= 1, "consistency SD"] = np.nan

# =========================
# 8. [新增] 计算 consensus Tiers
#    将 consensus score 四舍五入保留整数
# =========================
df["consensus Tiers"] = df["consensus score"].round(0)

# =========================
# 9. 计算 final rank
# 综合排序规则：
#   1) consensus Tiers 降序 (优先)
#   2) support count   降序
#   3) consistency SD  升序 (方差越小越好)
#   4) consensus score 降序 (终极平局打破者，保证稳定性)
# =========================
sort_df = df.copy()

# 为了避免 NaN 影响排序，建立临时列填充 NaN：
sort_df["_tier_tmp"] = sort_df["consensus Tiers"].fillna(-1)
sort_df["_support_tmp"] = sort_df["support count"].fillna(0)
sort_df["_sd_tmp"] = sort_df["consistency SD"].fillna(999999) # 缺失时设为极大值排在后面
sort_df["_score_tmp"] = sort_df["consensus score"].fillna(-1)

# 执行多级排序
sort_df = sort_df.sort_values(
    by=["_tier_tmp", "_support_tmp", "_sd_tmp", "_score_tmp"],
    ascending=[False, False, True, False]
).reset_index(drop=True)

# 生成最终排名 1, 2, 3...
sort_df["final rank"] = np.arange(1, len(sort_df) + 1)

# 把 rank 合并回原表 (通过 index 匹配后重新排序)
df = df.loc[sort_df.index.copy()].copy()
df["final rank"] = sort_df["final rank"].values

# 按照最终排名整理表格行顺序
df = df.sort_values("final rank").reset_index(drop=True)

# =========================
# 10. 保存结果
# =========================
# 输出1：推荐文件
df.to_excel(output_file1, index=False)

# 输出2：同内容另存一份
df.to_excel(output_file2, index=False)

print("✅ 计算与排名完成！")
print(f"输入文件: {input_file}")
print("\n【综合排名策略 (Tie-breaking)】:")
print("1. 优先按照 consensus Tiers (取整后分层) 降序排列")
print("2. Tiers 相同时，按照 support count 降序排列")
print("3. Support count 相同时，按照 consistency SD 升序排列 (越小越一致)")
print(f"\n结果已保存至:\n - {output_file1}\n - {output_file2}")