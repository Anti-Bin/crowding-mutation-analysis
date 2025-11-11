#!/usr/bin/env bash
set +e

# ============================================================
# Variant Processing and Annotation Pipeline
# ============================================================
# Author: [Bin Yu]
# Date: 2025-10-20
# Description:
#   This pipeline processes multiple VCF files to:
#     - extract flanking 3bp sequence around variants
#     - remove homopolymers
#     - decompose and normalize variants
#     - identify total and unique variants
#     - count and classify SNP/INDEL types
#     - annotate genomic features and coding impacts
# ============================================================

# --------------------------- CONFIG --------------------------
GENOME="/path/to/WBcel235.fa"
ORIGINAL_VCF_DIR="/path/to/original-vcf-files"
ANNOTATION_DIR="/path/to/original-annotation-files"
SCRIPT="/path/to/perl_code"
THREADS=4

# ------------------------ PREPARATION -------------------------
mkdir -p 3bp_3homodel \
    3bp_3homodel/vcf_decom \
    3bp_3homodel/vcf_decom/uniq \
    3bp_3homodel/vcf_decom/uniq/homodel_uniq_format \
    3bp_3homodel/vcf_decom/uniq/homodel_uniq_format/merge \
    3bp_3homodel/vcf_decom/uniq/homodel_uniq_format/specific \
    3bp_3homodel/vcf_decom/uniq/homodel_uniq_format/specific/count \

# -------------------- Extract 3bp sequences ----------------
echo "Extracting ±3bp FASTA sequences around variants..."
for id in *.vcf; do
    [[ -e "$id" ]] || { echo "No VCF files found!"; break; }
done
ls *.vcf | while read id; do
    perl -F"\t" -lane '
        unless(/^#/){
            my $x=$F[1]-4; my $y=$F[1]+3;
            my $a="$F[0]" . ":" . "$x" . "-" . "$y";
            print $F[0],"\t",$x,"\t",$y,"\t",$F[1],"\t",$a
        }' "$id" > 3bp_3homodel/${id%%.*}.3.bed

    bedtools getfasta -fi "$GENOME" \
                      -bed 3bp_3homodel/${id%%.*}.3.bed \
                      -fo 3bp_3homodel/${id%%.*}.3.fa.bed -tab
done

# --------------- Remove 3+ homopolymer sequences ------------
echo "Filtering variants from 3+ homopolymer sequences..."
cd 3bp_3homodel
for id in *.3.fa.bed; do
    base=${id%%.*}
    awk '$2 !~ /AAA|TTT|CCC|GGG/' "$id" > ${base}_3bp_3homodel.bed
    awk 'FNR==NR {a[$1]=$2; next} ($5 in a) {print}' \
        ${base}_3bp_3homodel.bed ${base}.3.bed \
        > ${base}.index.bed

    grep "^#" "$ORIGINAL_VCF_DIR/${base}.vcf" > ${base}.3bp.3homodel.vcf
    awk 'FNR==NR {a[$1_$4]; next} ($1_$2 in a) {print}' \
        ${base}.index.bed "$ORIGINAL_VCF_DIR/${base}.vcf" \
        >> ${base}.3bp.3homodel.vcf
done

# ---------------- Variant decomposition ---------------------
echo "Decomposing and normalizing variants..."
ls *.vcf | while read id; do
    vt decompose -s -o vcf_decom/${id%%.*}.decompose.vcf "$id"
    vt normalize -r "$GENOME" \
        -o vcf_decom/${id%%.*}.decompose_normalize.vcf \
        vcf_decom/${id%%.*}.decompose.vcf
done

# ---------------- Unique variant filtering ------------------
echo "Removing variants reported in corresponding F1..."
cd vcf_decom
ls cpr-4*.decompose_normalize.vcf | while read id; do
    base=${id%%.*}
    grep "^#" "$id" > uniq/${base}.3.uniq.homodel.vcf
    awk 'FNR==NR {a[$1_$2_$4_$5]; next} !($1_$2_$4_$5 in a) {print}' \
        cpr-4_G1.decompose_normalize.vcf "$id" | grep -v "^#" \
        >> uniq/${base}.3.uniq.homodel.vcf
done

ls N2*.decompose_normalize.vcf | while read id; do
    base=${id%%.*}
    grep "^#" "$id" > uniq/${base}.3.uniq.homodel.vcf
    awk 'FNR==NR {a[$1_$2_$4_$5]; next} !($1_$2_$4_$5 in a) {print}' \
        N2_G1.decompose_normalize.vcf "$id" | grep -v "^#" \
        >> uniq/${base}.3.uniq.homodel.vcf
done

# ---------------- Export query information ------------------
echo "Extracting query fields using bcftools..."
cd uniq
ls *.vcf | while read id; do
    base=${id%%.*}
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AF\t%DP\t%AO\t%TYPE\n' "$id" \
        > homodel_uniq_format/1.temp

    perl -lane 'unless(/^#/){$F[7]=~s/(.*)TYPE=(\w*)(.*)/$2/; print uc($F[7])}' "$id" \
        > homodel_uniq_format/2.temp

    paste -d"\t" homodel_uniq_format/1.temp homodel_uniq_format/2.temp \
        > homodel_uniq_format/${base}.homodel_uniq_formated

    rm homodel_uniq_format/1.temp homodel_uniq_format/2.temp
done

# ---------------- Merge and find common variants ------------
echo "Merging and identifying common variants..."
cd homodel_uniq_format
cat cpr-4_C* | sort | uniq > merge/cpr-4_10C.uniq
cat cpr-4_N* | sort | uniq > merge/cpr-4_10N.uniq
cat N2_N* | sort | uniq > merge/N2_10N.uniq
cat N2_C* | sort | uniq > merge/N2_10C.uniq

awk 'FNR==NR {a[$1_$2_$3_$4]; next} $1_$2_$3_$4 in a {print}' \
    merge/N2_10N.uniq merge/N2_10C.uniq > merge/N2_10C_10N.common

awk 'FNR==NR {a[$1_$2_$3_$4]; next} $1_$2_$3_$4 in a {print}' \
    merge/cpr-4_10N.uniq merge/cpr-4_10C.uniq > merge/cpr-4_10C_10N.common

# ------------ Remove variants from contrary group ------------
echo "Identifying group-specific variants..."
ls cpr-4_*_10* | while read id; do
    base=${id%%.*}
    awk 'FNR==NR {a[$1_$2_$3_$4]; next} !($1_$2_$3_$4 in a) {print}' \
        merge/cpr-4_10C_10N.common "$id" \
        > specific/${base}.homodel_uniq_specific
done

ls *N2_*_10* | while read id; do
    base=${id%%.*}
    awk 'FNR==NR {a[$1_$2_$3_$4]; next} !($1_$2_$3_$4 in a) {print}' \
        merge/N2_10C_10N.common "$id" \
        > specific/${base}.homodel_uniq_specific
done

# ---------------- Count variant quality metrics -------------
echo "Counting filtered variants..."
cd specific
ls *_specific | while read id; do
    base=${id%%.*}
    echo -ne "${base}\t" >> count/maf0.5_Q25_DP30_filter
    awk '{if ($6>=0.5 && $5>=25 && $7<=30 && $7>=5) {print}}' "$id" | wc -l \
        >> count/maf0.5_Q25_DP30_filter
done

# ---------------- Count variant types -----------------------
echo "Counting variant types..."
ls *_specific | while read id; do
    base=${id%%.*}
    for j in SNP MNP DEL INS; do
        echo -ne "${base}\t${j}\t" >> count/type_count_0.5_Q25_DP30_homo
        awk -v a="${j}" '{if ($6>=0.5 && $5>=25 && $7<=30 && $7>=5 && $10==a) {print}}' "$id" | wc -l \
            >> count/type_count_0.5_Q25_DP30_homo
    done
done

# ---------------- SNP subclassification --------------------
echo "Counting six SNP substitution types..."
for id in *_specific; do
    base=${id%%.*}
    awk '{if ($10=="SNP" && $6>=0.5 && $5>=25 && $7<=30 && $7>=5) {print $0}}' "$id" \
        > ${base}.snp_0.5_Q25_DP30_homo
done

for id in *.snp_0.5_Q25_DP30_homo; do
    base=${id%%.*}
    awk -v id="${base}" '
    BEGIN {
        types["T>C"]=0; types["T>G"]=0; types["T>A"]=0;
        types["C>G"]=0; types["C>T"]=0; types["C>A"]=0
    }
    {
        ref=toupper($3); alt=toupper($4);
        if ((ref=="A"&&alt=="G")||(ref=="T"&&alt=="C")) types["T>C"]++
        else if ((ref=="A"&&alt=="C")||(ref=="T"&&alt=="G")) types["T>G"]++
        else if ((ref=="A"&&alt=="T")||(ref=="T"&&alt=="A")) types["T>A"]++
        else if ((ref=="C"&&alt=="G")||(ref=="G"&&alt=="C")) types["C>G"]++
        else if ((ref=="C"&&alt=="T")||(ref=="G"&&alt=="A")) types["C>T"]++
        else if ((ref=="C"&&alt=="A")||(ref=="G"&&alt=="T")) types["C>A"]++
    }
    END {
        for (t in types) print id "\t" t "\t" types[t]
    }' "$id" >> count/SNP_change_count_0.5_Q25_DP30_homo
done

