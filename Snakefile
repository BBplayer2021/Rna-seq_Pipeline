# 导入必要的包
import os
from pathlib import Path
import pandas as pd

# 加载配置文件
configfile: "config.yaml"

# 读取metadata文件
metadata = pd.read_csv(config["metadata"], sep="\t")
SAMPLES = metadata["sample_id"].tolist()
RESULTS_DIR = config["file_params"]["results_dir"]

def is_paired_end(sample_id):
    """判断样本是否为双端测序"""
    sample_info = metadata[metadata["sample_id"] == sample_id].iloc[0]
    return sample_info["type"] == "paired"

def get_fastq(wildcards):
    """根据样本类型返回fastq文件"""
    sample_info = metadata[metadata["sample_id"] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        r1_pattern = config["naming_pattern"]["paired_end"]["r1"]
        r2_pattern = config["naming_pattern"]["paired_end"]["r2"]
        return {
            "r1": f"raw_data/{r1_pattern.format(sample=wildcards.sample)}{config['file_params']['fastq_ext']}", 
            "r2": f"raw_data/{r2_pattern.format(sample=wildcards.sample)}{config['file_params']['fastq_ext']}"
        }
    else:
        r1_pattern = config["naming_pattern"]["single_end"]["r1"]
        return {
            "r1": f"raw_data/{r1_pattern.format(sample=wildcards.sample)}{config['file_params']['fastq_ext']}"
        }

# 获取样本类型列表
PAIRED_SAMPLES = metadata[metadata["type"] == "paired"]["sample_id"].tolist()
SINGLE_SAMPLES = metadata[metadata["type"] == "single"]["sample_id"].tolist()

# 定义最终输出文件
rule all:
    input:
        f"{RESULTS_DIR}/multiqc/multiqc_report.html",
        expand(f"{RESULTS_DIR}/counts/{{sample}}_counts.txt", sample=SAMPLES),
        f"{RESULTS_DIR}/qc/fastqc_report.html",
        f"{RESULTS_DIR}/star/alignment_stats.txt"

# 导入各个步骤的规则
include: "workflow/quality_control.rules"
include: "workflow/alignment.rules"
include: "workflow/quantification.rules"
include: "workflow/report.rules" 