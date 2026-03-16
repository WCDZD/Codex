#!/usr/bin/env bash
set -euo pipefail

# ==============================
# 参数
# ==============================
if [ $# -lt 5 ]; then
  echo "Usage: $0 <samples.csv> <ref_root> <output_root> <threads_per_sample> <memory_per_sample>"
  exit 1
fi

SAMPLES_CSV="$1"
REF_ROOT="$2"
OUTPUT_ROOT="$3"
THREADS_PER="$4"       # 每个样本用多少线程
MEM_PER="$5"           # 每个样本申请多少内存

# ==============================
# 自动计算最大并行
# ==============================
TOTAL_MEM_GB=$(awk '/MemTotal/ {printf("%d\n", $2 / 1024 / 1024 + 0.5)}' /proc/meminfo)           # 总内存 (GB)
AVAIL_MEM_GB=$(awk '/MemAvailable/ {printf("%d\n", $2 / 1024 / 1024 + 0.5)}' /proc/meminfo)       # 当前可用内存 (GB)
TOTAL_CORES=$(nproc --all)

# 每个样本需要的内存（转成 GB 数字）
MEM_GB=${MEM_PER%G}   # 去掉 G，得到数字，如 80G -> 80

# 计算最大并行样本数（取内存和CPU限制的最小值）
MAX_JOBS_MEM=$(( AVAIL_MEM_GB / MEM_GB ))                  # 内存允许的最多线程
MAX_JOBS_CPU=$(( TOTAL_CORES / THREADS_PER ))               # CPU 允许的最多线程
MAX_JOBS=$(( MAX_JOBS_MEM < MAX_JOBS_CPU ? MAX_JOBS_MEM : MAX_JOBS_CPU ))

# 至少跑 1 个，最大不超过样本数
MAX_JOBS=${MAX_JOBS:-1}
if [ "$MAX_JOBS" -lt 1 ]; then MAX_JOBS=1; fi

echo "系统总内存: ${TOTAL_MEM_GB}GB，可用: ${AVAIL_MEM_GB}GB"
echo "总核心: ${TOTAL_CORES}"
echo "每个样本: ${THREADS_PER} 线程, ${MEM_PER} 内存"
echo "自动计算最大并发样本数: $MAX_JOBS"
echo "样本表路径: $SAMPLES_CSV"
echo "样本表内容:"
cat "$SAMPLES_CSV"


# 创建临时目录
TMP_CSV_DIR="${OUTPUT_ROOT}/tmp_csv"
mkdir -p "$TMP_CSV_DIR"

# 提取数据行到临时文件（每行一个临时 CSV）
HEADER="SampleID,tumor_fq_1_path,tumor_fq_2_path,normal_fq_1_path,normal_fq_2_path"

# 生成临时单行 CSV（带表头）
awk -v header="$HEADER" -F, 'NR==1 {next}  # 跳过原表头
NR>1 {gsub(/[\r]/,""); 
  sample=$1
  tmp_csv = "'"$TMP_CSV_DIR"'/" sample ".csv"
  print header > tmp_csv
  print $0 >> tmp_csv
  print tmp_csv
}' "$SAMPLES_CSV" > "${OUTPUT_ROOT}/tmp_csv_list.txt"


# 并行运行（用 --arg-file 传临时 CSV 路径） 加-u实时输出Log,但多样本同时跑会很乱
parallel -j "$MAX_JOBS" \
  --verbose \
  --joblog "${OUTPUT_ROOT}/parallel_joblog.txt" \
  bash /data/xuzl/WESpipeline/pipeline_docker.sh \
    {} \
    "$REF_ROOT" \
    "${OUTPUT_ROOT}/result" \
    "$THREADS_PER" \
    "$MEM_PER" \
  :::: "${OUTPUT_ROOT}/tmp_csv_list.txt"


# -------------------------- 批处理完成后的智能升级 --------------------------
echo "第一轮完成，开始判断是否构建 PoN 并升级 CNV/HRD"

# 计算有效样本数（跳过表头和空行）
SAMPLE_COUNT=$(awk -F, 'NR>1 && $1!="" && $1!~/#/' "$SAMPLES_CSV" | wc -l)
PON_CNV_MARK="${OUTPUT_ROOT}/00.PON_CNV_COMPLETE"

echo "检测到有效样本数: $SAMPLE_COUNT"
if [ "$SAMPLE_COUNT" -le 2 ]; then
  echo "样本数 <3，不构建 PoN，流程结束（已使用 hybrid mode 完成分析）"
elif [ -f "$PON_CNV_MARK" ]; then
  echo "检测到 PoN CNV已完成标记: $PON_CNV_MARK"
  echo "跳过 PoN 构建和第二轮重跑（已使用 PoN + allele-specific + 完整 HRD）"
else
  echo "样本数 ≥3，开始构建 CNVkit PoN 并升级 CNVkit + HRD"
  # 1. 构建 PoN 这里根据这个批次构建一个pon放在总目录，同时流程也会自动识别REF_ROOT下是否有cnvkit_pon文件用作全局通用Pon（如果以后有更大人群的通用pon文件可以直接用的话比用NormalBam自己构建好)
  echo "构建 PoN reference 和 PureCN normalDB..."
  PON_DIR="${OUTPUT_ROOT}/cnvkit_pon"
  PURECN_NORMALDB_DIR="${OUTPUT_ROOT}/purecn_normalDB"
  mkdir -p "$PON_DIR" "$PURECN_NORMALDB_DIR"

  if [ ! -f "${PON_DIR}/reference.cnn" ]; then
  bash /data/xuzl/WESpipeline/build_cnvkit_pon.sh \
    "$SAMPLES_CSV" \
    "$REF_ROOT" \
    "${OUTPUT_ROOT}/result" \
    "$PON_DIR" \
    "$THREADS_PER"
  else
    echo "${PON_DIR}/reference.cnn 已存在，跳过"
  fi

  # 2. 构建 PureCN normalDB 和 mapping_bias
  INTERVAL_FILE="${REF_ROOT}/purecn/intervals_hg38.txt"
  NORMAL_BAM_LIST="${PURECN_NORMALDB_DIR}/normal_bams.list" 
  find "${OUTPUT_ROOT}/result" -name "*_normal.sorted.dedup.recal.bam" > "$NORMAL_BAM_LIST"

  echo "找到 $(wc -l < "$NORMAL_BAM_LIST") 个 normal BAM，list 文件：$NORMAL_BAM_LIST"
  head "$NORMAL_BAM_LIST"

# 创建 normal coverage 输出目录
  NORMAL_COV_DIR="${PURECN_NORMALDB_DIR}/normal_coverage"
  mkdir -p "$NORMAL_COV_DIR"

# 获取 Coverage.R 路径
  COVERAGE_SCRIPT=$(/opt/conda/bin/mamba run -n purecn Rscript -e "cat(system.file('extdata', 'Coverage.R', package='PureCN'))")

# 运行 Coverage.R（关键：--bam 后面跟 .list 文件）
  /opt/conda/bin/mamba run -n purecn Rscript "$COVERAGE_SCRIPT" \
    --out-dir "$NORMAL_COV_DIR" \
    --bam "$NORMAL_BAM_LIST" \
    --intervals "$INTERVAL_FILE" \
    --cores 8

  LOESS_LIST="${PURECN_NORMALDB_DIR}/normal_loess.list"  

  find "$NORMAL_COV_DIR" -name "*_loess.txt.gz" > "$LOESS_LIST"

  echo "loess 文件列表（$LOESS_LIST）："
  cat "$LOESS_LIST"

# 确保文件可读
  chmod 644 "$LOESS_LIST"
  NORMALDB_SCRIPT=$(/opt/conda/bin/mamba run -n purecn Rscript -e "cat(system.file('extdata', 'NormalDB.R', package='PureCN'))")

# 构建 normalDB
  /opt/conda/bin/mamba run -n purecn Rscript "$NORMALDB_SCRIPT" \
    --out-dir "$PURECN_NORMALDB_DIR" \
    --coverage-files "$LOESS_LIST" \
    --genome hg38 \
    --assay agilent_v8


  echo "normalDB 构建结果："
  ls -l "$PURECN_NORMALDB_DIR"/*_hg38.rds
  PURECN_READY_MARK="${OUTPUT_ROOT}/00.PURECN_READY"
  touch "$PURECN_READY_MARK"
  echo "创建 PureCN 就绪标记: $PURECN_READY_MARK"
  echo "CNVkit PoN 和 PureCN normalDB 构建完成"


  # 2. 删除所有样本的 CNVkit 和 HRD .SUCCESS（强制重跑）
  echo "删除 CNVkit 和 HRD 标记文件，准备重跑..."
  find "${OUTPUT_ROOT}/result" -path "*/08.CNV/.SUCCESS" -delete
  find "${OUTPUT_ROOT}/result" -path "*/09.TMB_MSI_HRD/.SUCCESS" -delete 2>/dev/null || true

  # 4. 创建升级完成标记
  touch "$PON_CNV_MARK"
  echo "已使用PON更新CNV结果，创建标记: $PON_CNV_MARK"
  

  # 3. 第二轮并行重跑（只重跑 CNVkit + HRD + 后续）
  echo "第二轮并行重跑：使用 PoN 升级 CNVkit 和完整 HRD"
  parallel -j "$MAX_JOBS" \
    --verbose \
    --joblog "${OUTPUT_ROOT}/parallel_joblog_round2.txt" \
    bash /data/xuzl/WESpipeline/pipeline_docker.sh \
      {} \
      "$REF_ROOT" \
      "${OUTPUT_ROOT}/result" \
      "$THREADS_PER" \
      "$MEM_PER" \
    :::: "${OUTPUT_ROOT}/tmp_csv_list.txt"

  
fi

# -------------------------- 最终结果汇总 --------------------------
echo "开始最终结果汇总..."
bash /data/xuzl/WESpipeline/summarize_results.sh "${OUTPUT_ROOT}/result"

# 清理临时文件
rm -rf "$TMP_CSV_DIR" "${OUTPUT_ROOT}/tmp_csv_list.txt"

echo "======================================================================"
echo "所有样本批量处理完成！结果在 $OUTPUT_ROOT 下"
echo "======================================================================"