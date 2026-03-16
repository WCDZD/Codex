#!/usr/bin/env bash
set -euo pipefail

# ==============================
# 参数解析
# ==============================
if [ $# -lt 5 ]; then
  echo "Usage: $0 <samples.csv> <ref_root> <output_root> <threads> <memory>"
  echo "Note: samples.csv must be comma-separated (CSV format)"
  echo "Example: $0 /data/samples.tsv /data/reference /data/results 24 \"80G\""
  exit 1
fi

SAMPLES_TSV="$1"
REF_ROOT="$2"          # 需要传入的参考根目录
RESULT_ROOT="$3"
THREADS="$4"
MEM="$5"

# ==============================
# 所有参考路径基于 REF_ROOT 
# ==============================
REF="${REF_ROOT}/GRCh38/Homo_sapiens_assembly38.fasta"
KNOWN="${REF_ROOT}/knownsites"
ANNOVAR_DB="${REF_ROOT}/annovar/humandb"
VEP_CACHE="${REF_ROOT}/vep_cache" 
BED="${REF_ROOT}/WES_bed/Agilent_SureSelect_V8_GRCh38_Regions_clean.bed"
REFFLAT="${REF_ROOT}/WES_bed/refFlat.txt"
ACCESS_CHR="${REF_ROOT}/WES_bed/access.hg38.bed"
HLAHD_DIR="/opt/hlahd"   # 镜像内固定



MSI_SITES="${REF_ROOT}/msisensor/hg38.microsatellite.list"
HRD_SCORE_SITES="${REF_ROOT}/hrd/SNP6.hg38.interval_list"
CHROMOSOMES="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"



# ==============================
# 辅助函数：检查并执行步骤（带跳过机制）
# ==============================
run_step() {
  local step_dir="$1"      # 步骤输出目录
  local step_name="$2"     # 步骤名称（用于日志显示）
  local success_file="${step_dir}/.SUCCESS"

  if [ -f "$success_file" ]; then
    echo "======================================================================"
    echo "[SKIP] $step_name 已完成（检测到 $success_file），跳过"
    echo "======================================================================"
    return 0  # 表示已完成，可跳过
  else
    echo "======================================================================"
    echo "[RUN] 正在执行: $step_name"
    echo "======================================================================"
    return 1  # 需要运行
  fi
}

mark_step_done() {
  local step_dir="$1"
  local step_name="$2"
  local success_file="${step_dir}/.SUCCESS"

  touch "$success_file"
  echo "======================================================================"
  echo "[DONE] $step_name 完成，生成标记文件: $success_file"
  echo "======================================================================"
}


# ==============================
# 读取样本表并循环处理
# ==============================
IFS=$'\t'
while IFS=, read -r SAMPLE TUMOR_R1 TUMOR_R2 NORMAL_R1 NORMAL_R2; do
  # 去除可能的尾部 \r（Windows 换行符）
  SAMPLE=${SAMPLE%$'\r'}
  TUMOR_R1=${TUMOR_R1%$'\r'}
  TUMOR_R2=${TUMOR_R2%$'\r'}
  NORMAL_R1=${NORMAL_R1%$'\r'}
  NORMAL_R2=${NORMAL_R2%$'\r'}

  # 跳过空行、注释行，以及表头行（常见表头关键词）
  [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# || "$SAMPLE" =~ ^SampleID|^sample_id|^SAMPLE ]] && continue

  echo "======================================================================"
  echo "Processing Sample: $SAMPLE"
  echo "Tumor  FASTQ: $TUMOR_R1  $TUMOR_R2"
  echo "Normal FASTQ: $NORMAL_R1  $NORMAL_R2"
  echo "======================================================================"

  RESULT_PATH="${RESULT_ROOT}/${SAMPLE}"
  TRIMMED_DIR="${RESULT_PATH}/02.Trimmed"
  TMP_DIR="${RESULT_PATH}/tmp"

  mkdir -p "${RESULT_PATH}/01.QC/raw_fastqc" \
           "${RESULT_PATH}/01.QC/multiqc" \
           "${RESULT_PATH}/02.Trimmed" \
           "${RESULT_PATH}/03.Alignment" \
           "${RESULT_PATH}/04.Mutect2" \
           "${RESULT_PATH}/05.Annotation" \
           "${RESULT_PATH}/06.HLA" \
           "${RESULT_PATH}/07.Neoantigen" \
           "${RESULT_PATH}/08.CNV" \
           "${RESULT_PATH}/09.TMB_MSI_HRD" \
           "${RESULT_PATH}/logs" \
           "${TMP_DIR}"


  exec > >(tee -a "${RESULT_PATH}/logs/${SAMPLE}_pipeline.log") 2>&1
  
  FASTQ_TMP="${TMP_DIR}/fastq"
  mkdir -p "$FASTQ_TMP"

  export TMPDIR="${TMP_DIR}"
  export MPLCONFIGDIR="${TMP_DIR}/matplotlib"
  export XDG_CACHE_HOME="${TMP_DIR}/xdg_cache"
  export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=${TMP_DIR}/java"
  mkdir -p "${TMP_DIR}/matplotlib" \
         "${TMP_DIR}/xdg_cache" \
         "${TMP_DIR}/java" \
         "${TMP_DIR}/vep" \
         "${TMP_DIR}/pvacseq"
  
  # 防止fontconfig报错
  export FONTCONFIG_PATH="${TMP_DIR}/fontconfig"
  export FONTCONFIG_FILE="${TMP_DIR}/fontconfig/fonts.conf"
  mkdir -p "${TMP_DIR}/fontconfig"
  mkdir -p "${TMP_DIR}/fontconfig"
  touch "${TMP_DIR}/fontconfig/fonts.conf"

  # -------------------------- 1. QC --------------------------
  if run_step "${RESULT_PATH}/01.QC" "01. QC (FastQC + MultiQC)"; then
    true  # 跳过
  else
    fastqc -t 12 $TUMOR_R1 $TUMOR_R2 $NORMAL_R1 $NORMAL_R2 \
      -o "${RESULT_PATH}/01.QC/raw_fastqc" --nogroup

    multiqc "${RESULT_PATH}/01.QC/raw_fastqc" \
      -o "${RESULT_PATH}/01.QC/multiqc" \
     -n "${SAMPLE}_multiqc_report.html" \
      --title "WES QC - ${SAMPLE}"

    mark_step_done "${RESULT_PATH}/01.QC" "01. QC"
  fi

  # -------------------------- 2. Trim --------------------------
  #get_basename() {
  #  local name=$(basename "$1")
  #  name="${name%.gz}"; name="${name%.fastq}"; name="${name%.fq}"
  #  echo "$name"
  #}

  #NORMAL_BASE=$(get_basename "$NORMAL_R1")
  #TUMOR_BASE=$(get_basename "$TUMOR_R1")
  NORMAL_BASE="${SAMPLE}_normal"
  TUMOR_BASE="${SAMPLE}_tumor"
  if run_step "$TRIMMED_DIR" "02. Trim (fastp)"; then
    # 即使跳过，也要确保更新后的路径正确（因为后面依赖）
    NORMAL_R1="${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R1.fq.gz"
    NORMAL_R2="${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R2.fq.gz"
    TUMOR_R1="${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R1.fq.gz"
    TUMOR_R2="${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R2.fq.gz"
  else
    fastp -i "$NORMAL_R1" -I "$NORMAL_R2" \
      -o "${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R1.fq.gz" \
      -O "${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R2.fq.gz" \
      --detect_adapter_for_pe --qualified_quality_phred 20 \
      --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
      --length_required 36 --thread "$THREADS" \
      --html "${TRIMMED_DIR}/normal_fastp.html" \
      --json "${TRIMMED_DIR}/normal_fastp.json"

    fastp -i "$TUMOR_R1" -I "$TUMOR_R2" \
      -o "${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R1.fq.gz" \
      -O "${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R2.fq.gz" \
      --detect_adapter_for_pe --qualified_quality_phred 20 \
      --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
      --length_required 36 --thread "$THREADS" \
      --html "${TRIMMED_DIR}/tumor_fastp.html" \
      --json "${TRIMMED_DIR}/tumor_fastp.json"

    # 更新路径为 trimmed
    NORMAL_R1="${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R1.fq.gz"
    NORMAL_R2="${TRIMMED_DIR}/${NORMAL_BASE}_trimmed_R2.fq.gz"
    TUMOR_R1="${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R1.fq.gz"
    TUMOR_R2="${TRIMMED_DIR}/${TUMOR_BASE}_trimmed_R2.fq.gz"

    mark_step_done "$TRIMMED_DIR" "02. Trim"
  fi

  # -------------------------- 3. Alignment --------------------------
  if run_step "${RESULT_PATH}/03.Alignment" "03. Alignment + BQSR"; then
    true
  else
    bwa mem -t "$THREADS" -M -R "@RG\tID:${SAMPLE}_tumor\tSM:${SAMPLE}_tumor\tPL:ILLUMINA\tLB:WES" \
    "$REF" "$TUMOR_R1" "$TUMOR_R2" | \
    samtools view -@ $((THREADS/4)) -Sb - > "${RESULT_PATH}/03.Alignment/${SAMPLE}_tumor.unsorted.bam"
    bwa mem -t "$THREADS" -M -R "@RG\tID:${SAMPLE}_normal\tSM:${SAMPLE}_normal\tPL:ILLUMINA\tLB:WES" \
    "$REF" "$NORMAL_R1" "$NORMAL_R2" | \
    samtools view -@ $((THREADS/4)) -Sb - > "${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.unsorted.bam"
    for TYPE in tumor normal; do
      gatk --java-options "-Xmx${MEM}" SortSam \
      -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.unsorted.bam" \
      -O "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.bam" \
      --SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR "$TMP_DIR"
      gatk --java-options "-Xmx${MEM}" MarkDuplicates \
      -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.bam" \
      -O "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.dedup.bam" \
      -M "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.metrics.txt" \
      --CREATE_INDEX true --TMP_DIR "$TMP_DIR"
    done
  # -------------------------- 3.5. BQSR --------------------------
    for TYPE in tumor normal; do
      gatk --java-options "-Xmx${MEM}" BaseRecalibrator \
      -R "$REF" -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.dedup.bam" \
      --known-sites "${KNOWN}/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
      --known-sites "${KNOWN}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
      --known-sites "${KNOWN}/af-only-gnomad.hg38.vcf.gz" \
      -O "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.recal.table"
    
      gatk --java-options "-Xmx${MEM}" ApplyBQSR \
      -R "$REF" -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.dedup.bam" \
      --bqsr-recal-file "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.recal.table" \
      -O "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.dedup.recal.bam"

      samtools index "${RESULT_PATH}/03.Alignment/${SAMPLE}_${TYPE}.sorted.dedup.recal.bam"
    done

    mark_step_done "${RESULT_PATH}/03.Alignment" "03&04. Alignment + BQSR"
  fi


  

  # -------------------------- 4. Mutect2 --------------------------
  if run_step "${RESULT_PATH}/04.Mutect2" "04. Mutect2"; then
    true
  else
    gatk --java-options "-Xmx${MEM}" Mutect2 \
    -R "$REF" \
    -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_tumor.sorted.dedup.recal.bam" \
    -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.sorted.dedup.recal.bam" \
    -normal "${SAMPLE}_normal" \
    --panel-of-normals "${KNOWN}/somatic-hg38_1000g_pon.hg38.vcf.gz" \
    --germline-resource "${KNOWN}/af-only-gnomad.hg38.vcf.gz" \
    -O "${RESULT_PATH}/04.Mutect2/${SAMPLE}.unfiltered.vcf.gz" \
    --tmp-dir "$TMP_DIR"
    
    # ========== 4.2 GetPileupSummaries（tumor 和 normal）==========
    TUMOR_PILEUP="${RESULT_PATH}/04.Mutect2/${SAMPLE}_tumor.pileups.table"
    NORMAL_PILEUP="${RESULT_PATH}/04.Mutect2/${SAMPLE}_normal.pileups.table"

    echo "[Mutect2] Generating pileup summaries..."
    gatk --java-options "-Xmx${MEM}" GetPileupSummaries \
      -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_tumor.sorted.dedup.recal.bam" \
      -V "${KNOWN}/af-only-gnomad.hg38.vcf.gz" \
      -L "${BED}" \
      -R "$REF" \
      -O "$TUMOR_PILEUP" \
      --tmp-dir "$TMP_DIR"

    gatk --java-options "-Xmx${MEM}" GetPileupSummaries \
      -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.sorted.dedup.recal.bam" \
      -V "${KNOWN}/af-only-gnomad.hg38.vcf.gz" \
      -L "${BED}" \
      -R "$REF" \
      -O "$NORMAL_PILEUP" \
      --tmp-dir "$TMP_DIR"

    # ========== 4.3 CalculateContamination ==========
    CONTAM_TABLE="${RESULT_PATH}/04.Mutect2/${SAMPLE}.contamination.table"
    SEGMENTS_TABLE="${RESULT_PATH}/04.Mutect2/${SAMPLE}.segments.table"

    echo "[Mutect2] Calculating contamination..."
    gatk --java-options "-Xmx${MEM}" CalculateContamination \
      -I "$TUMOR_PILEUP" \
      -matched "$NORMAL_PILEUP" \
      -O "$CONTAM_TABLE" \
      --tumor-segmentation "$SEGMENTS_TABLE" \
      --tmp-dir "$TMP_DIR"

    # ========== 4.4 FilterMutectCalls==========
    gatk FilterMutectCalls -R "$REF" \
    -V "${RESULT_PATH}/04.Mutect2/${SAMPLE}.unfiltered.vcf.gz" \
    -O "${RESULT_PATH}/04.Mutect2/${SAMPLE}.filtered.vcf.gz" \
    --contamination-table "$CONTAM_TABLE" \
    --tumor-segmentation "$SEGMENTS_TABLE"

    bcftools view -f PASS "${RESULT_PATH}/04.Mutect2/${SAMPLE}.filtered.vcf.gz" -Oz \
    -o "${RESULT_PATH}/04.Mutect2/${SAMPLE}.PASS.vcf.gz"
    tabix -p vcf "${RESULT_PATH}/04.Mutect2/${SAMPLE}.PASS.vcf.gz"

    mark_step_done "${RESULT_PATH}/04.Mutect2" "04. Mutect2"
  fi

  # ==============================
  # 4.5 Germline Variant Calling (胚系突变调用)
  # ==============================
  GERMLINE_DIR="${RESULT_PATH}/04.Germline"
  mkdir -p "$GERMLINE_DIR"

  if run_step "$GERMLINE_DIR" "04.5 Germline Variant Calling (HaplotypeCaller)"; then
    echo "Germline calling 已完成，跳过"
  else
    echo "Running GATK HaplotypeCaller for germline variants from normal sample..."

    gatk --java-options "-Xmx${MEM}" HaplotypeCaller \
      -R "$REF" \
      -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.sorted.dedup.recal.bam" \
      -O "${GERMLINE_DIR}/${SAMPLE}_germline.raw.vcf.gz" \
      --native-pair-hmm-threads "$THREADS" \
      --standard-min-confidence-threshold-for-calling 20 \
      --output-mode EMIT_VARIANTS_ONLY \
      --tmp-dir "$TMP_DIR"

    # Hard-filter germline variants (GATK best practices for germline)
    gatk --java-options "-Xmx${MEM}" VariantFiltration \
      -V "${GERMLINE_DIR}/${SAMPLE}_germline.raw.vcf.gz" \
      -O "${GERMLINE_DIR}/${SAMPLE}_germline.filtered.vcf.gz" \
      --filter-name "LowQual" \
      --filter-expression "QUAL < 30.0" \
      --filter-name "LowDepth" \
      --filter-expression "DP < 10" \
      --filter-name "StrandBias" \
      --filter-expression "FS > 60.0" \
      --filter-name "LowMappingQual" \
      --filter-expression "MQ < 40.0" \
      --filter-name "ReadPosEndDist" \
      --filter-expression "ReadPosEndDist > 3.0" \
      --genotype-filter-name "LowGQ" \
      --genotype-filter-expression "GQ < 20" \
      --tmp-dir "$TMP_DIR"

    # Extract PASS variants
    bcftools view -f PASS "${GERMLINE_DIR}/${SAMPLE}_germline.filtered.vcf.gz" -Oz \
      -o "${GERMLINE_DIR}/${SAMPLE}_germline.PASS.vcf.gz"

    tabix -p vcf "${GERMLINE_DIR}/${SAMPLE}_germline.PASS.vcf.gz"

    # 最终胚系 VCF 输出
    FINAL_GERMLINE_VCF="${GERMLINE_DIR}/${SAMPLE}_germline.PASS.vcf.gz"
    echo "胚系突变 VCF 已生成：$FINAL_GERMLINE_VCF"
    echo "   - 原始调用：${GERMLINE_DIR}/${SAMPLE}_germline.raw.vcf.gz"
    echo "   - 过滤后 PASS：$FINAL_GERMLINE_VCF"

    mark_step_done "$GERMLINE_DIR" "04.5 Germline Variant Calling"
  fi
  # ==============================
  # 4.6 Germline Variant Annotation (ANNOVAR)
  # ==============================
  GERMLINE_ANNO_DIR="${GERMLINE_DIR}/anno"
  mkdir -p "$GERMLINE_ANNO_DIR"
  if run_step "$GERMLINE_ANNO_DIR" "04.6 Germline Variant Annotation (ANNOVAR)"; then
    true
  else
    echo "Annotating germline variants with Annovar..."

    GERMLINE_PASS_VCF="${GERMLINE_DIR}/${SAMPLE}_germline.PASS.vcf.gz"

    table_annovar.pl "${GERMLINE_PASS_VCF}" "$ANNOVAR_DB" \
      -vcfinput \
      -buildver hg38 \
      -out "${GERMLINE_ANNO_DIR}/${SAMPLE}_germline_anno" \
      -remove \
      -protocol refGene,clinvar_20240917,gnomad312_genome,dbnsfp47a,avsnp147,exac03 \
      -operation g,f,f,f,f,f \
      -nastring . \
      -polish -otherinfo
    FINAL_GERMLINE_ANNO="${GERMLINE_ANNO_DIR}/${SAMPLE}_germline_anno.hg38_multianno.csv"

    # 添加样本列
    awk -F'\t' -v sample="$SAMPLE" '
      BEGIN {OFS=","}
      NR==1 {print "Sample", $0}
      NR>1  {print sample, $0}
    ' "${GERMLINE_ANNO_DIR}/${SAMPLE}_germline_anno.hg38_multianno.txt" > "$FINAL_GERMLINE_ANNO"

    echo "Germline ANNOVAR 注释完成：$FINAL_GERMLINE_ANNO"
    mark_step_done "$GERMLINE_ANNO_DIR" "04.6 Germline Variant Annotation (ANNOVAR)"
  fi

  # -------------------------- 5. Annotation --------------------------
  if run_step "${RESULT_PATH}/05.Annotation" "05. Annotation (Annovar + VEP)"; then
    true
  else
    convert2annovar.pl -format vcf4 -allsample --outfile "${RESULT_PATH}/05.Annotation/${SAMPLE}" \
    "${RESULT_PATH}/04.Mutect2/${SAMPLE}.PASS.vcf.gz"
    TUMOR_AVINPUT=$(ls "${RESULT_PATH}/05.Annotation/${SAMPLE}"*tumor.avinput | head -1)

    table_annovar.pl "$TUMOR_AVINPUT" "$ANNOVAR_DB" \
    -buildver hg38 -out "${RESULT_PATH}/05.Annotation/${SAMPLE}" -remove \
    -protocol refGene,clinvar_20240917,gnomad312_genome,dbnsfp47a,avsnp147,exac03 \
    -operation g,f,f,f,f,f \
    -nastring . -polish -otherinfo

    ANNOVAR_TXT="${RESULT_PATH}/05.Annotation/${SAMPLE}.hg38_multianno.txt"
    ANNOVAR_FINAL="${RESULT_PATH}/05.Annotation/${SAMPLE}_final.hg38_multianno.csv"

    if [ ! -f "$ANNOVAR_TXT" ]; then
      echo "ERROR: ANNOVAR output txt not found: $ANNOVAR_TXT"
      exit 1
    fi

    # 把 txt 转成 csv 并添加 Tumor_Sample_Barcode 列
    awk -F'\t' -v sample="$SAMPLE" '
      BEGIN {OFS=","}
      NR==1 {print "Tumor_Sample_Barcode", $0}
      NR>1  {print sample, $0}
    ' "$ANNOVAR_TXT" > "$ANNOVAR_FINAL"


    # VEP for pVACseq
    VEP_VCF="${RESULT_PATH}/05.Annotation/${SAMPLE}.vep.vcf"
    VEP_VCF_GZ="${VEP_VCF}.gz"

    vep --input_file "${RESULT_PATH}/04.Mutect2/${SAMPLE}.PASS.vcf.gz" \
      --output_file "$VEP_VCF" --format vcf --vcf \
      --symbol --hgvs --protein --biotype --canonical --tsl --mane_select \
      --fasta "$REF" --offline --cache --dir_cache "$VEP_CACHE" \
      --dir_plugins "${VEP_CACHE}/Plugins" \
      --plugin Frameshift --plugin Wildtype \
      --fork "$THREADS" --force_overwrite
    
    bgzip -c "$VEP_VCF" > "$VEP_VCF_GZ.tmp"  # 临时文件

    # 用 bcftools sort 强制排序
    bcftools sort "$VEP_VCF_GZ.tmp" -Oz -o "$VEP_VCF_GZ"
    tabix -p vcf "$VEP_VCF_GZ"
    # 清理临时
    rm -f "$VEP_VCF_GZ.tmp"


    mark_step_done "${RESULT_PATH}/05.Annotation" "05. Annotation"
  fi

  # -------------------------- 6. HLA Typing --------------------------
  if run_step "${RESULT_PATH}/06.HLA" "06. HLA Typing"; then
    true
  else
    process_fastq() {
      local input="$1"
      case "$input" in
      *.gz) gunzip -c "$input" > "${FASTQ_TMP}/$(basename "$input" .gz)" && echo "${FASTQ_TMP}/$(basename "$input" .gz)" ;;
      *) echo "$input" ;;
      esac
    }
    NORMAL_R1_P=$(process_fastq "$NORMAL_R1")
    NORMAL_R2_P=$(process_fastq "$NORMAL_R2")
    TUMOR_R1_P=$(process_fastq "$TUMOR_R1")
    TUMOR_R2_P=$(process_fastq "$TUMOR_R2")

    # OptiType (Class I)
    OptiTypePipeline.py -i $NORMAL_R1_P $NORMAL_R2_P --dna -v -o "${RESULT_PATH}/06.HLA" -p "${SAMPLE}_opti_normal"
    OptiTypePipeline.py -i $TUMOR_R1_P $TUMOR_R2_P --dna -v -o "${RESULT_PATH}/06.HLA" -p "${SAMPLE}_opti_tumor"

  # HLA-HD (Class I+II)
    ulimit -n 4096
    mkdir -p "${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_normal" "${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_tumor"
    hlahd.sh -t "$THREADS" -m 50 -f "${HLAHD_DIR}/freq_data" \
      "$NORMAL_R1_P" "$NORMAL_R2_P" \
      "${HLAHD_DIR}/HLA_gene.split.txt" "${HLAHD_DIR}/dictionary" \
      "${SAMPLE}_normal" "${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_normal"
  

    hlahd.sh -t "$THREADS" -m 50 -f "${HLAHD_DIR}/freq_data" \
      "$TUMOR_R1_P" "$TUMOR_R2_P" \
      "${HLAHD_DIR}/HLA_gene.split.txt" "${HLAHD_DIR}/dictionary" \
      "${SAMPLE}_tumor" "${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_tumor"
    
    mark_step_done "${RESULT_PATH}/06.HLA" "06. HLA Typing"
  fi

  # ---- HLA 结果提取 & LOH 判断（保持原逻辑）----
  
  # ========== HLA结果提取、对比、LOH判断（最终版：三字段转两字段） ==========
  echo "Extracting HLA alleles from HLA-HD and converting to pVACseq format..."
  # OptiType Class I (带 HLA- 前缀，两字段)
  OPT_NORMAL_CLASSI=$(awk 'NR==2 {print "HLA-"$2",HLA-"$3",HLA-"$4",HLA-"$5",HLA-"$6",HLA-"$7}' ${RESULT_PATH}/06.HLA/${SAMPLE}_opti_normal_result.tsv)
  OPT_TUMOR_CLASSI=$(awk 'NR==2 {print "HLA-"$2",HLA-"$3",HLA-"$4",HLA-"$5",HLA-"$6",HLA-"$7}' ${RESULT_PATH}/06.HLA/${SAMPLE}_opti_tumor_result.tsv)

  echo "OptiType Class I Normal: $OPT_NORMAL_CLASSI"
  echo "OptiType Class I Tumor : $OPT_TUMOR_CLASSI"

  # HLA-HD Class II (DRB1 无前缀 + DQ/DP 配对无前缀)
  extract_class_ii() {
    local file=$1
    local raw=$(grep -v "Not typed" "$file" | \
      awk '{for(i=3;i<=NF;i++) if($i ~ /^HLA-[A-Z]/ && $i !~ /-$/) print $i}')

    local two=$(echo "$raw" | sed 's/\([A-Z*]*:[0-9]*\):[0-9]*/\1/g')

    local drb1=$(echo "$two" | grep 'HLA-DRB1' | sed 's/HLA-//g' | tr '\n' ',' | sed 's/,$//')

    local dqa=$(echo "$two" | grep 'HLA-DQA1' | sed 's/HLA-//g' | sort)
    local dqb=$(echo "$two" | grep 'HLA-DQB1' | sed 's/HLA-//g' | sort)
    local dq_pairs=""
    if [ -n "$dqa" ] && [ -n "$dqb" ]; then
      dq_pairs=$(paste <(echo "$dqa") <(echo "$dqb") | awk '{print $1"-"$2}' | tr '\n' ',' | sed 's/,$//')
    fi

    local dpa=$(echo "$two" | grep 'HLA-DPA1' | sed 's/HLA-//g' | sort)
    local dpb=$(echo "$two" | grep 'HLA-DPB1' | sed 's/HLA-//g' | sort)
    local dp_pairs=""
    if [ -n "$dpa" ] && [ -n "$dpb" ]; then
      dp_pairs=$(paste <(echo "$dpa") <(echo "$dpb") | awk '{print $1"-"$2}' | tr '\n' ',' | sed 's/,$//')
    fi

    local all=""
    [ -n "$drb1" ] && all="$drb1"
    [ -n "$dq_pairs" ] && all="$all,$dq_pairs"
    [ -n "$dp_pairs" ] && all="$all,$dp_pairs"

    echo "$all" | sed 's/^,//; s/,$//; s/ //g'
  }

  HLahd_NORMAL_RESULT=${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_normal/${SAMPLE}_normal/result/${SAMPLE}_normal_final.result.txt
  HLahd_TUMOR_RESULT=${RESULT_PATH}/06.HLA/${SAMPLE}_hlahd_tumor/${SAMPLE}_tumor/result/${SAMPLE}_tumor_final.result.txt

  HLahd_NORMAL_CLASSII=$(extract_class_ii "$HLahd_NORMAL_RESULT")
  HLahd_TUMOR_CLASSII=$(extract_class_ii "$HLahd_TUMOR_RESULT")

  echo "HLA-HD Class II Normal: $HLahd_NORMAL_CLASSII"
  echo "HLA-HD Class II Tumor : $HLahd_TUMOR_CLASSII"

  # LOH 判断（分别判断 Class I 和 Class II）
  use_class_i=""
  if [ "$OPT_NORMAL_CLASSI" = "$OPT_TUMOR_CLASSI" ]; then
    echo "No Class I LOH detected → Using Normal Class I"
    use_class_i="$OPT_NORMAL_CLASSI"
  else
    echo "Potential Class I LOH detected → Using Tumor Class I"
    use_class_i="$OPT_TUMOR_CLASSI"
  fi

  use_class_ii=""
  if [ "$HLahd_NORMAL_CLASSII" = "$HLahd_TUMOR_CLASSII" ]; then
    echo "No Class II LOH detected → Using Normal Class II"
    use_class_ii="$HLahd_NORMAL_CLASSII"
  else
    echo "Potential Class II LOH detected → Using Tumor Class II"
    use_class_ii="$HLahd_TUMOR_CLASSII"
  fi

  # 组合最终 alleles
  FINAL_HLA_ALLELES="$use_class_i"
  [ -n "$use_class_ii" ] && FINAL_HLA_ALLELES="$FINAL_HLA_ALLELES,$use_class_ii"
  FINAL_HLA_ALLELES=$(echo "$FINAL_HLA_ALLELES" | sed 's/^,//; s/,$//; s/ //g')
  echo "Final HLA alleles for pVACseq (OptiType Class I + HLA-HD Class II): $FINAL_HLA_ALLELES"


  # -------------------------- 7. pVACseq --------------------------
  if run_step "${RESULT_PATH}/07.Neoantigen/${SAMPLE}" "07. pVACseq"; then
    true
  else
    rm -rf "${RESULT_PATH}/07.Neoantigen/${SAMPLE}"
    mkdir -p "${RESULT_PATH}/07.Neoantigen/${SAMPLE}"
    chmod 777 "${RESULT_PATH}/07.Neoantigen/${SAMPLE}"

    TEMP_HOME="/tmp/pv_${SAMPLE}"
    rm -rf "$TEMP_HOME/.local/share/mhcflurry"
    mkdir -p "$TEMP_HOME/.local/share/mhcflurry"
    ln -sf /opt/mhcflurry_data "$TEMP_HOME/.local/share/mhcflurry/4"
    #
    #MHCnuggetsI 
    env -i PATH="/usr/local/bin:/usr/bin:/bin" \
           HOME="$TEMP_HOME" \
           OMP_NUM_THREADS=1 \
           TF_NUM_INTEROP_THREADS=1 \
           TF_NUM_INTRAOP_THREADS=1 \
           MKL_NUM_THREADS=1 \
           /usr/local/bin/pvacseq run \
      "${RESULT_PATH}/05.Annotation/${SAMPLE}.vep.vcf.gz" \
      "${SAMPLE}_tumor" \
      "$FINAL_HLA_ALLELES" \
      NetMHCpanEL MHCflurry BigMHC_EL NetMHCIIpanEL\
      "${RESULT_PATH}/07.Neoantigen/${SAMPLE}" \
      -e1 8,9,10,11 -e2 13,15,17,19 \
      -t 1 \
      --binding-threshold 500 \
      --percentile-threshold 2 \
      --top-score-metric lowest \
      --keep-tmp-files \
      --iedb-install-directory /opt/iedb

    # 清理
    rm -rf "$TEMP_HOME"
    
    mark_step_done "${RESULT_PATH}/07.Neoantigen/${SAMPLE}" "07. pVACseq"
  fi
  

  

  # -------------------------- 8. CNVkit --------------------------
  TUMOR_BAM="${RESULT_PATH}/03.Alignment/${SAMPLE}_tumor.sorted.dedup.recal.bam"
  NORMAL_BAM="${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.sorted.dedup.recal.bam"
  CNV_DIR="${RESULT_PATH}/08.CNV"

  # -------------------------- Germline SNP for BAF (optional but recommended for allele-specific CNV) --------------------------
  GERMLINE_VCF="${CNV_DIR}/germline/${SAMPLE}_germline.vcf.gz"
  mkdir -p "${CNV_DIR}/germline"
  echo "test"
  if run_step "${CNV_DIR}/germline" "Germline SNP calling"; then
    true
  else
    gatk --java-options "-Xmx${MEM}" HaplotypeCaller \
      -R "$REF" \
      -I "$NORMAL_BAM" \
      -O "$GERMLINE_VCF" \
      --standard-min-confidence-threshold-for-calling 20 \
      --output-mode EMIT_VARIANTS_ONLY \
      --sample-name "${SAMPLE}_normal"

    # 过滤 heterozygous SNPs（可选，更精确）
    bcftools view -v snps -g het "$GERMLINE_VCF" -Oz -o "${GERMLINE_VCF%.vcf.gz}.het.vcf.gz"
    tabix -p vcf "${GERMLINE_VCF%.vcf.gz}.het.vcf.gz"

    mark_step_done "${CNV_DIR}/germline" "Germline SNP calling"
  fi

  if run_step "$CNV_DIR" "08. CNVkit"; then
    true
  else
    GLOBAL_PON_REF="${REF_ROOT}/cnvkit_pon/reference.cnn"  # 全局的PON
    NORMALIZED_ROOT="${RESULT_ROOT%/}"          # 先去掉结尾斜杠
    BATCH_ROOT="${NORMALIZED_ROOT%/*}"          # 再取上级目录
    BATCH_PON_REF="${BATCH_ROOT}/cnvkit_pon/reference.cnn"

    if [ -f "$BATCH_PON_REF" ]; then
      PON_REF="$BATCH_PON_REF"
      echo "[CNV] Found batch-specific PoN: $PON_REF (highest priority)"
    elif [ -f "$GLOBAL_PON_REF" ]; then
      PON_REF="$GLOBAL_PON_REF"
      echo "[CNV] Found global PoN: $PON_REF"
    else
      PON_REF=""
      echo "[CNV] No PoN found (batch or global), fallback to hybrid mode"
    fi

    if [ -n "$PON_REF" ]; then
      echo "[CNV] Found PoN reference: $PON_REF"
      echo "[CNV] Running CNVkit with PoN (allele-specific capable)"
      cnvkit.py batch "$TUMOR_BAM" \
        --reference "$PON_REF" \
        --output-dir "$CNV_DIR" \
        -p "$THREADS" \
        --scatter --diagram

    else
      echo "[CNV] No PoN reference found at $PON_REF"
      echo "[CNV] Fallback to hybrid mode with matched normal"
      cnvkit.py batch "$TUMOR_BAM" --normal "$NORMAL_BAM" \
        --targets "$BED" --fasta "$REF" --access "$ACCESS_CHR" \
        --annotate "$REFFLAT" --output-reference "${CNV_DIR}/${SAMPLE}_reference.cnn" \
        --output-dir "$CNV_DIR" -p "$THREADS" --scatter --diagram --method hybrid
    fi

    cnvkit.py genemetrics "${CNV_DIR}/${SAMPLE}_tumor.sorted.dedup.cnr" \
      -s "${CNV_DIR}/${SAMPLE}_tumor.sorted.dedup.cns" \
      -o "${CNV_DIR}/${SAMPLE}_gene_metrics.txt"

    cnvkit.py call "${CNV_DIR}/${SAMPLE}_tumor.sorted.dedup.cnr" \
      -o "${CNV_DIR}/${SAMPLE}_tumor.call.cns"

    cnvkit.py genemetrics "${CNV_DIR}/${SAMPLE}_tumor.sorted.dedup.cnr" \
      -s "${CNV_DIR}/${SAMPLE}_tumor.call.cns" \
      --threshold 0 \
      -o "${CNV_DIR}/${SAMPLE}_gene_absolute_cn.txt"
  
    cnvkit.py call "${CNV_DIR}/${SAMPLE}_tumor.call.cns" \
      -v "${GERMLINE_VCF%.vcf.gz}.het.vcf.gz" \
      -o "${CNV_DIR}/${SAMPLE}_tumor.allele_specific.call.cns"

    cnvkit.py genemetrics "${CNV_DIR}/${SAMPLE}_tumor.sorted.dedup.cnr" \
      -s "${CNV_DIR}/${SAMPLE}_tumor.allele_specific.call.cns" \
      --threshold 0 \
      -o "${CNV_DIR}/${SAMPLE}_gene_major_minor_absolute_cn.txt"


    NORMAL_CNR="${CNV_DIR}/${SAMPLE}_normal.sorted.dedup.cnr"
    if [ ! -f "$NORMAL_CNR" ]; then
      echo "[CNVkit] 为 normal bam 生成独立 .cnr 文件（用于后续 PureCN normalDB）"
      cnvkit.py coverage "$NORMAL_BAM" "$BED" \
        -o "$NORMAL_CNR" \
        -p "$THREADS"
    else
      echo "[CNVkit] normal .cnr 已存在，跳过"
    fi


    mark_step_done "$CNV_DIR" "08. CNVkit"
  fi

  # -------------------------- 8.5 Purity & Ploidy (PureCN) - 仅当 batch PoN 存在时运行 --------------------------
  PURITY_DIR="${CNV_DIR}/PurityPloidy"
  PURECN_OUT="${PURITY_DIR}/${SAMPLE}"
  mkdir -p "$PURECN_OUT"

  INTERVAL_FILE="${REF_ROOT}/purecn/intervals_hg38.txt"

  # 重新计算 batch normalDB 路径（与批量脚本一致）
  BATCH_ROOT="${RESULT_ROOT%/}"          # 去掉结尾斜杠
  BATCH_ROOT="${BATCH_ROOT%/*}"          # 取上级目录
  PURECN_READY_MARK="${BATCH_ROOT}/00.PURECN_READY"

  if [ ! -f "$PURECN_READY_MARK" ]; then
  # 第一轮或 normalDB 未准备好：跳过 PureCN，输出 NA
  echo "[PureCN] PureCN ready mark not found ($PURECN_READY_MARK), skipping PureCN (first round or normalDB not built)."
  cat > "${PURITY_DIR}/${SAMPLE}_purity_ploidy.txt" <<EOF
Sample	Purity	Ploidy	Sex	Status
$SAMPLE	NA	NA	NA	Skipped (waiting for batch normalDB)
EOF
    # 不标记 .SUCCESS
  else
  # 第二轮：normalDB 已就绪，执行完整 PureCN
    if run_step "$PURITY_DIR" "08.5 Purity & Ploidy (PureCN)"; then
      echo "[PureCN] Already completed, skipping."
    else
      echo "[PureCN] PureCN ready mark found, starting PureCN analysis..."

    # 此时可以安全查找 normalDB（标记存在说明已构建）
    BATCH_NORMALDB=$(find "$BATCH_ROOT/purecn_normalDB" -name "normalDB_*_hg38.rds" | head -1)
      # ========== 1. 为 PureCN 生成专用 VCF（包含 somatic + germline het SNPs）==========
      PURECN_VCF="${PURECN_OUT}/${SAMPLE}_purecn_input.vcf.gz"
      PURECN_VCF_UNFILTERED="${PURECN_OUT}/${SAMPLE}_purecn_unfiltered.vcf.gz"

      if [ ! -f "$PURECN_VCF" ]; then
        echo "[PureCN] Running dedicated Mutect2 for PureCN (with germline genotyping)..."

        gatk --java-options "-Xmx${MEM}" Mutect2 \
          -R "$REF" \
          -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_tumor.sorted.dedup.recal.bam" \
          -I "${RESULT_PATH}/03.Alignment/${SAMPLE}_normal.sorted.dedup.recal.bam" \
          -normal "${SAMPLE}_normal" \
          --panel-of-normals "${KNOWN}/somatic-hg38_1000g_pon.hg38.vcf.gz" \
          --germline-resource "${KNOWN}/af-only-gnomad.hg38.vcf.gz" \
          --genotype-germline-sites true \
          --genotype-pon-sites true \
          --interval-padding 75 \
          -O "$PURECN_VCF_UNFILTERED" \
          --tmp-dir "$TMP_DIR"

      echo "[PureCN] Filtering dedicated VCF (more lenient for germline SNPs)..."
      gatk FilterMutectCalls \
        -R "$REF" \
        -V "$PURECN_VCF_UNFILTERED" \
        -O "$PURECN_VCF" 

      tabix -p vcf -f "$PURECN_VCF"
      echo "[PureCN] Dedicated VCF ready: $PURECN_VCF"
    else
      echo "[PureCN] Using existing dedicated VCF: $PURECN_VCF"
    fi

    # ========== 2. 计算 tumor coverage loess（与 normalDB 完全匹配）==========
    TUMOR_LOESS_DIR="${PURECN_OUT}/coverage"
    mkdir -p "$TUMOR_LOESS_DIR"

    TUMOR_LOESS=$(find "$TUMOR_LOESS_DIR" -name "*_loess.txt.gz" | head -1)

    if [ ! -f "$TUMOR_LOESS" ]; then
      echo "[PureCN] Calculating tumor coverage loess..."
      COVERAGE_SCRIPT=$(/opt/conda/bin/mamba run -n purecn Rscript -e "cat(system.file('extdata', 'Coverage.R', package='PureCN'))")

      /opt/conda/bin/mamba run -n purecn Rscript "$COVERAGE_SCRIPT" \
        --out-dir "$TUMOR_LOESS_DIR" \
        --bam "$TUMOR_BAM" \
        --intervals "$INTERVAL_FILE" \
        --cores "$THREADS" \
        --force

      TUMOR_LOESS=$(find "$TUMOR_LOESS_DIR" -name "*_loess.txt.gz" | head -1)
      echo "[PureCN] Tumor loess file: $TUMOR_LOESS"
    else
      echo "[PureCN] Using existing tumor loess: $TUMOR_LOESS"
    fi

    # ========== 3. 运行 PureCN ==========
    echo "[PureCN] Running PureCN.R..."
    PURECN_SCRIPT=$(/opt/conda/bin/mamba run -n purecn Rscript -e "cat(system.file('extdata', 'PureCN.R', package='PureCN'))")

    /opt/conda/bin/mamba run -n purecn Rscript "$PURECN_SCRIPT" \
      --out "$PURECN_OUT" \
      --tumor "$TUMOR_LOESS" \
      --sampleid "$SAMPLE" \
      --vcf "$PURECN_VCF" \
      --normaldb "$BATCH_NORMALDB" \
      --intervals "$INTERVAL_FILE" \
      --fun-segmentation PSCBS \
      --genome hg38 \
      --force \
      --seed 123 \
      --min-base-quality 20 \
      --min-supporting-reads 3 \
      --error 0.002 \
      --max-segments 500

    # ========== 4. 提取 purity/ploidy 结果 ==========
    if [ -f "${PURECN_OUT}/${SAMPLE}.csv" ]; then
      # PureCN 主输出文件通常是 ${SAMPLE}.csv
      RAW_PURITY=$(awk -F',' 'NR==2 {gsub(/"/,"",$2); print $2}' "${PURECN_OUT}/${SAMPLE}.csv")
      RAW_PLOIDY=$(awk -F',' 'NR==2 {gsub(/"/,"",$3); print $3}' "${PURECN_OUT}/${SAMPLE}.csv")
      SEX=$(awk -F',' 'NR==2 {gsub(/"/,"",$NF); print $NF}' "${PURECN_OUT}/${SAMPLE}.csv" | cut -d' ' -f1)
      STATUS="Success"
      PURITY=$(printf "%.2f" "$RAW_PURITY")
      PLOIDY=$(printf "%.3f" "$RAW_PLOIDY")
      echo "[PureCN] Success: Purity=$PURITY, Ploidy=$PLOIDY, Sex=$SEX"
    else
      PURITY="NA"; PLOIDY="NA"; SEX="NA"; STATUS="Failed"
      echo "[PureCN] Failed: No output CSV found."
    fi

    cat > "${PURITY_DIR}/${SAMPLE}_purity_ploidy.txt" <<EOF
Sample	Purity	Ploidy	Sex	Status
$SAMPLE	$PURITY	$PLOIDY	$SEX	$STATUS
EOF

    mark_step_done "$PURITY_DIR" "08.5 Purity & Ploidy (PureCN)"
  fi
fi



  # -------------------------- 9. TMB + MSI + HRD --------------------------
  TMB_MSI_HRD_DIR="${RESULT_PATH}/09.TMB_MSI_HRD"
  MSI_DIR="${TMB_MSI_HRD_DIR}/MSI"
  TMB_DIR="${TMB_MSI_HRD_DIR}/TMB"
  HRD_DIR="${TMB_MSI_HRD_DIR}/HRD"


  if run_step "$TMB_MSI_HRD_DIR" "09. TMB + MSI + HRD"; then
    true
  else
    mkdir -p "$MSI_DIR" "$TMB_DIR" "$HRD_DIR"

    PURECN_SUMMARY="${PURITY_DIR}/${SAMPLE}_purity_ploidy.txt"
    if [ -f "$PURECN_SUMMARY" ] && grep -q "Success" "$PURECN_SUMMARY"; then
      PURECN_PURITY=$(awk 'NR==2 {print $2}' "$PURECN_SUMMARY")
      PURECN_PLOIDY=$(awk 'NR==2 {print $3}' "$PURECN_SUMMARY")
      PURECN_SEX=$(awk 'NR==2 {print $4}' "$PURECN_SUMMARY")
      PURECN_STATUS=$(awk 'NR==2 {print $5}' "$PURECN_SUMMARY")
      echo "[HRD] Found successful PureCN result: Purity=$PURECN_PURITY, Ploidy=$PURECN_PLOIDY, Sex=$PURECN_SEX"
      HRD_NOTE="Full HRD score using allele-specific segments + PureCN ploidy"
    else
      PURECN_PURITY="NA"; PURECN_PLOIDY="NA"; PURECN_SEX="NA"
      echo "[HRD] No successful PureCN result found, falling back to original HRD calculation."
      HRD_NOTE="Full HRD score using allele-specific segments (no PureCN ploidy)"
    fi

    # ========== 9.1 MSI (msisensor-pro) ==========
    if run_step "$MSI_DIR" "MSI"; then
      true
    else
      echo "[09.1] Running msisensor-pro for MSI detection..."
      msisensor-pro msi \
        -d "${REF_ROOT}/msisensor/hg38.microsatellite.list" \
        -n "$NORMAL_BAM" \
        -t "$TUMOR_BAM" \
        -o "${MSI_DIR}/${SAMPLE}_msi" \
        -b "$THREADS" \
        -c 20 \
        -z 1
      mark_step_done "$MSI_DIR" "MSI"
    fi

    # 提取结果（第二行最后一列是 % somatic）
    MSI_SCORE=$(awk 'NR==2 {print $(NF)}' "${MSI_DIR}/${SAMPLE}_msi" 2>/dev/null || echo "NA")
    if awk "BEGIN {exit !($MSI_SCORE>10)}"; then
      MSI_STATUS="MSI-H"
    else
      MSI_STATUS="MSS"
    fi
    [ -z "$MSI_SCORE" ] && MSI_SCORE="NA" && MSI_STATUS="NA"

    # ========== 9.2 TMB (Tumor Mutational Burden) ==========
    echo "[09.2] Calculating TMB..."
    SOMATIC_VCF="${RESULT_PATH}/04.Mutect2/${SAMPLE}.PASS.vcf.gz"
    NEO_FINAL="${RESULT_PATH}/07.Neoantigen/${SAMPLE}/MHC_Class_I/${SAMPLE}_tumor.MHC_I.filtered.tsv"  # Class I 过滤结果
    NEO_FINAL_II="${RESULT_PATH}/07.Neoantigen/${SAMPLE}/MHC_Class_II/${SAMPLE}_tumor.MHC_II.filtered.tsv"  # Class II 过滤结果（如果存在）

    # 计算捕获区域大小（Mb）
    CAPTURE_SIZE_MB=$(awk '{sum += $3-$2} END {printf "%.2f", sum/1000000}' "$BED")
    # ========== TMB ==========
    if [ -f "$SOMATIC_VCF" ] && [ "$CAPTURE_SIZE_MB" != "0.00" ]; then
      TOTAL_MUT=$(bcftools view -f PASS "$SOMATIC_VCF" | grep -v "^#" | wc -l)
      TMB=$(awk "BEGIN {printf \"%.2f\", $TOTAL_MUT / $CAPTURE_SIZE_MB}")
      # TMB 状态判断
      if awk "BEGIN {exit !($TMB >= 10)}"; then
        TMB_STATUS="TMB-High"
      else
        TMB_STATUS="TMB-Low"
      fi
    else
      TOTAL_MUT="NA"
      TMB="NA"
      TMB_STATUS="NA"
    fi
    # ========== TNB (Tumor Neoantigen Burden) ==========
    if [ -f "$NEO_FINAL" ]; then
      # Class I neoantigen 数量（跳过表头）
      NEO_COUNT_I=$(awk 'NR>1' "$NEO_FINAL" | wc -l)
    else
      NEO_COUNT_I=0
    fi
    if [ -f "$NEO_FINAL_II" ]; then
      # Class II neoantigen 数量（跳过表头）
      NEO_COUNT_II=$(awk 'NR>1' "$NEO_FINAL_II" | wc -l)
    else
      NEO_COUNT_II=0
    fi

    NEO_COUNT_TOTAL=$((NEO_COUNT_I + NEO_COUNT_II))

    if [ "$(awk "BEGIN {print ($CAPTURE_SIZE_MB > 0)}")" = "1" ] && [ "$NEO_COUNT_TOTAL" -gt 0 ]; then
      TNB=$(awk "BEGIN {printf \"%.2f\", $NEO_COUNT_TOTAL / $CAPTURE_SIZE_MB}")
    else
      TNB="0.00"
    fi

    echo "[TNB] Class I neoantigens: $NEO_COUNT_I"
    echo "[TNB] Class II neoantigens: $NEO_COUNT_II"
    echo "[TNB] Total neoantigens: $NEO_COUNT_TOTAL"
    echo "[TNB] Tumor Neoantigen Burden (TNB): $TNB per Mb"


    # ========== 9.4 HRD (scarHRD with allele-specific) ==========
    echo "[09.4] Calculating HRD..."
    PON_REF="${REF_ROOT}/cnvkit_pon/reference.cnn"
    CNS_FILE="${CNV_DIR}/${SAMPLE}_tumor.allele_specific.call.cns"

    if [ -f "$CNS_FILE" ]; then
      if [ -f "$PON_REF" ]; then
        echo "[HRD] Detected PoN reference: $PON_REF"
        echo "[HRD] Using scarHRD for full HRD score (LOH + LST + TAI, allele-specific available)"

        echo "[HRD] Running scarHRD analysis..."
        /opt/conda/bin/mamba run -n hrd_env Rscript - <<EOF
library(scarHRD)
library(data.table)
cat("Processing sample:", "$SAMPLE", "\n")

cns <- fread("$CNS_FILE", data.table = FALSE)
cns_complete <- cns[!is.na(cns\$cn1) & !is.na(cns\$cn2), ]

prepare_scarhrd_data <- function(cnvkit_data, sample_id, ploidy_val) {
  seg_data <- data.frame(
    sample = sample_id,
    chromosome = gsub("chr", "", cnvkit_data\$chromosome),
    start = as.integer(cnvkit_data\$start),
    end = as.integer(cnvkit_data\$end),
    total_cn = as.integer(cnvkit_data\$cn),
    A_cn = as.integer(cnvkit_data\$cn1),
    B_cn = as.integer(cnvkit_data\$cn2),
    stringsAsFactors = FALSE
  )
  seg_data\$ploidy <- if (is.na(ploidy_val) || ploidy_val == "NA") median(seg_data\$total_cn, na.rm = TRUE) else as.numeric(ploidy_val)
  seg_data\$chromosome <- ifelse(seg_data\$chromosome %in% c("X", "x"), "23",
                                ifelse(seg_data\$chromosome %in% c("Y", "y"), "24", seg_data\$chromosome))
  standard_chr <- c(as.character(1:22), "23", "24")
  seg_data <- seg_data[seg_data\$chromosome %in% standard_chr, ]
  seg_data <- seg_data[order(as.numeric(seg_data\$chromosome), seg_data\$start), ]
  return(seg_data)
}

scarhrd_input <- prepare_scarhrd_data(cns_complete, "$SAMPLE", "$PURECN_PLOIDY")
tmp_file <- tempfile(fileext = ".seg")
write.table(scarhrd_input, tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

result <- scar_score(tmp_file, reference = "grch38", seqz = FALSE)
hrd_scores <- data.frame(
  Sample = "$SAMPLE",
  HRD_Score = as.numeric(result[1, "HRD-sum"]),
  HRD_LOH = as.numeric(result[1, "HRD"]),
  LST = as.numeric(result[1, "LST"]),
  Telomeric_AI = as.numeric(result[1, "Telomeric AI"]),
  stringsAsFactors = FALSE
)

output_file <- "${HRD_DIR}/${SAMPLE}_scarHRD.txt"
write.table(hrd_scores, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Results saved to:", output_file, "\n")
cat("Used ploidy:", scarhrd_input\$ploidy[1], "\n")
EOF
        # 提取结果
        RESULTS_FILE="${HRD_DIR}/${SAMPLE}_scarHRD.txt"
        if [[ -f "$RESULTS_FILE" ]]; then
            HRD_SCORE=$(awk 'NR==2 {print $2}' "$RESULTS_FILE")
            LOH_SCORE=$(awk 'NR==2 {print $3}' "$RESULTS_FILE")
            LST_SCORE=$(awk 'NR==2 {print $4}' "$RESULTS_FILE")
            TAI_SCORE=$(awk 'NR==2 {print $5}' "$RESULTS_FILE")
            echo "[HRD] HRD scores: Total=$HRD_SCORE, LOH=$LOH_SCORE, LST=$LST_SCORE, TAI=$TAI_SCORE"
            if [[ "$HRD_SCORE" != "NA" ]]; then
              if awk -v hrd="$HRD_SCORE" 'BEGIN {exit !(hrd >= 42)}'; then
                HRD_STATUS="HRD-Positive"
              else
                HRD_STATUS="HRD-Negative"
              fi
            else
              HRD_STATUS="Failed"
            fi
            
            HRD_NOTE="Full HRD score using allele-specific segments"
        else
            HRD_SCORE="NA"
            LOH_SCORE="NA"
            LST_SCORE="NA"
            TAI_SCORE="NA"
            HRD_STATUS="Failed"
            HRD_NOTE="scarHRD analysis failed"
        fi
      else
        echo "[HRD] No PoN reference found. Skipping scarHRD (TAI unreliable in single-sample hybrid mode)."
        echo "[HRD] Calculating simplified HRD (LOH + LST only) as fallback."
        # 简化版 awk HRD
        awk '
        BEGIN { loh = 0; lst = 0; prev_chr = ""; prev_end = 0; prev_cn = -1 }
        NR > 1 {
          chr = $1; gsub("chr","",chr); start = $2; end = $3; cn = $5; len = end - start + 1;
          if (cn < 2 && len > 15000000) loh++;
          if (prev_chr == chr && prev_cn != cn && (start - prev_end) > 10000000) lst++;
          prev_chr = chr; prev_end = end; prev_cn = cn;
        }
        END {
          hrd = loh + lst;
          print "Sample\tHRD_Score\tLOH\tLST\tTAI" > "'${HRD_DIR}/${SAMPLE}_simpleHRD.txt'"
          print "'$SAMPLE'\t" hrd "\t" loh "\t" lst "\tNA" >> "'${HRD_DIR}/${SAMPLE}_simpleHRD.txt'"
        }' "$CNS_FILE"

        HRD_SCORE=$(awk 'NR==2 {print $2}' "${HRD_DIR}/${SAMPLE}_simpleHRD.txt")
        LOH_SCORE=$(awk 'NR==2 {print $3}' "${HRD_DIR}/${SAMPLE}_simpleHRD.txt")
        LST_SCORE=$(awk 'NR==2 {print $4}' "${HRD_DIR}/${SAMPLE}_simpleHRD.txt")
        TAI_SCORE="NA"
        if awk "BEGIN {exit !($HRD_SCORE >= 42)}"; then
          HRD_STATUS="HRD-Positive"
        else
          HRD_STATUS="HRD-Negative"
        fi
        HRD_NOTE="Simplified HRD score (LOH + LST only). No PoN reference available, TAI not calculated."
      fi
    else
      echo "[HRD] CNV call file not found: $CNS_FILE. Skipping HRD calculation."
      HRD_SCORE="NA"; LOH_SCORE="NA"; LST_SCORE="NA"; TAI_SCORE="NA"; HRD_STATUS="NA"
      HRD_NOTE="HRD not calculated (CNV call file missing)."
    fi

    # ========== Summary (tabular format) ==========
    SUMMARY_FILE="${TMB_MSI_HRD_DIR}/summary.txt"

    # 直接写入
    cat > "$SUMMARY_FILE" <<-EOF
Sample	TMB	Mutation_Count	Capture_Size_Mb	TNB	Neoantigen_Count_I	Neoantigen_Count_II	Neoantigen_Count_Total	MSI_Score	MSI_Status	HRD_Score	LOH	LST	TAI	HRD_Status	HRD_Note
EOF

    cat >> "$SUMMARY_FILE" <<-EOF
$SAMPLE	$TMB	$TOTAL_MUT	$CAPTURE_SIZE_MB	$TNB	$NEO_COUNT_I	$NEO_COUNT_II	$NEO_COUNT_TOTAL	$MSI_SCORE	$MSI_STATUS	$HRD_SCORE	$LOH_SCORE	$LST_SCORE	$TAI_SCORE	$HRD_STATUS	$HRD_NOTE
EOF

    echo "09. TMB_MSI_HRD tabular summary generated: $SUMMARY_FILE"
    echo "   Ready for multi-sample merging (e.g., cat */09.TMB_MSI_HRD/summary.txt > all_samples.tsv)"

  
    mark_step_done "$TMB_MSI_HRD_DIR" "09. TMB + MSI + HRD"
  fi

  echo "Sample $SAMPLE 完成！"
  

done < "$SAMPLES_TSV"

echo "======================================================================"
echo "所有样本处理完毕！"
echo "======================================================================"