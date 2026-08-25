import pandas as pd

def main():
    # 读取数据
    df = pd.read_excel("merge_with_scores.xlsx")
    
    # 多条件排序：
    # consensus Tiers 降序 -> support count 降序 -> consistency SD 升序 -> consensus score 降序
    df_sorted = df.sort_values(
        by=["consensus Tiers", "support count", "consistency SD", "consensus score"],
        ascending=[False, False, True, False]
    )
    
    # 重置索引
    df_sorted = df_sorted.reset_index(drop=True)
    total = len(df_sorted)
    
    # 生成排名列（从1开始）
    df_sorted["Rank"] = range(1, total + 1)
    # 计算排名占总体的百分数（%）
    df_sorted["Rank_Percent"] = (df_sorted["Rank"] / total * 100).round(2)
    
    # 保存结果
    df_sorted.to_excel("merge_with_scores_ranked.xlsx", index=False)
    print(f"排名完成！共 {total} 条数据，结果已保存至 merge_with_scores_ranked.xlsx")
    print(df_sorted.head(10))

if __name__ == "__main__":
    main()