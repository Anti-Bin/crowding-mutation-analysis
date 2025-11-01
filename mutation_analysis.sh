#!/usr/bin/env bash
set -euo pipefail

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
#     - identify unique variants
#     - count and classify SNP/INDEL types
#     - annotate genomic features and coding impacts
# ============================================================

# --------------------------- CONFIG --------------------------
GENOME="WBcel235.fa"
ORIGINAL_VCF_DIR="/path/to/original"
ANNOTATION_DIR="/path/to/original/anno"
THREADS=4

# ------------------------ PREPARATION -------------------------
mkdir -p AT_del count 3bp_3homodel vcf_decom uniq \
         homodel_uniq_format merge common_C_NC \
         SNP_type SNP_count anno anno_count/anno_count/count

# -------------------- 1. Extract 3bp sequences ----------------
echo "Extracting ±3bp FASTA sequences around variants..."
ls *.vcf | while read id; do
    perl -F"\t" -lane '
        unless(/^#/){
            my $x=$F[1]-4; my $y=$F[1]+3;
            my $a="$F[0]" . ":" . "$x" . "-" . "$y";
            print $F[0],"\t",$x,"\t",$y,"\t",$F[1],"\t",$a
        }' "$id" > AT_del/${id%%.*}.3.bed

    bedtools getfasta -fi "$GENOME" \
                      -bed AT_del/${id%%.*}.3.bed \
                      -fo AT_del/${id%%.*}.3.fa.bed -tab
done

# --------------- 2. Remove 3+ homopolymer sequences ------------
echo "Filtering variants from 3+ homopolymer sequences..."
for id in *.3.fa.bed; do
    base=${id%%.*}
    echo -ne "${base}\t" >> count/3homo_del_3bp
    awk '$2 !~ /AAA|TTT|CCC|GGG/' "$id" | wc -l >> count/3homo_del_3bp

    awk '$2 !~ /AAA|TTT|CCC|GGG/' "$id" > 3bp_3homodel/${base}_3bp_3homodel.bed
    awk 'FNR==NR {a[$1]=$2; next} ($5 in a) {print}' \
        3bp_3homodel/${base}_3bp_3homodel.bed ${base}.3.bed \
        > 3bp_3homodel/${base}.index.bed

    grep "^#" "$ORIGINAL_VCF_DIR/${base}.vcf" > 3bp_3homodel/${base}.3bp.3homodel.vcf
    awk 'FNR==NR {a[$1_$4]; next} ($1_$2 in a) {print}' \
        3bp_3homodel/${base}.index.bed "$ORIGINAL_VCF_DIR/${base}.vcf" \
        >> 3bp_3homodel/${base}.3bp.3homodel.vcf
done

# ---------------- 3. Variant decomposition ---------------------
echo "Decomposing and normalizing variants..."
ls *.vcf | while read id; do
    vt decompose -s -o vcf_decom/${id%%.*}.decompose.vcf "$id"
    vt normalize -r "$GENOME" \
        -o vcf_decom/${id%%.*}.decompose_normalize.vcf \
        vcf_decom/${id%%.*}.decompose.vcf
done

# ---------------- 4. Unique variant filtering ------------------
echo "Removing variants reported in corresponding F1..."
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

# ---------------- 5. Export query information ------------------
echo "Extracting query fields using bcftools..."
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

# ---------------- 6. Merge and find common variants ------------
echo "Merging and identifying common variants..."
cat cpr-4_C* | sort | uniq > merge/cpr-4_10C.uniq
cat cpr-4_N* | sort | uniq > merge/cpr-4_10N.uniq
cat N2_N* | sort | uniq > merge/N2_10N.uniq
cat N2_C* | sort | uniq > merge/N2_10C.uniq

awk 'FNR==NR {a[$1_$2_$3_$4]; next} $1_$2_$3_$4 in a {print}' \
    N2_10N.uniq N2_10C.uniq > N2_10C_10N.common

awk 'FNR==NR {a[$1_$2_$3_$4]; next} $1_$2_$3_$4 in a {print}' \
    cpr-4_10N.uniq cpr-4_10C.uniq > cpr-4_10C_10N.common

# ------------ 7. Remove variants from contrary group ------------
echo "Identifying group-specific variants..."
ls cpr-4_*_10-* | while read id; do
    base=${id%%.*}
    awk 'FNR==NR {a[$1_$2_$3_$4]; next} !($1_$2_$3_$4 in a) {print}' \
        cpr-4_10C_10N.common "$id" \
        > common_C_NC/${base}.homodel_uniq_specific
done

ls *N2_*_10-* | while read id; do
    base=${id%%.*}
    awk 'FNR==NR {a[$1_$2_$3_$4]; next} !($1_$2_$3_$4 in a) {print}' \
        N2_10C_10N.common "$id" \
        > common_C_NC/${base}.homodel_uniq_specific
done

# ---------------- 8. Count variant quality metrics -------------
echo "Counting filtered variants..."
ls *_specific | while read id; do
    base=${id%%.*}
    echo -ne "${base}\t" >> count/maf0.2_Q25_DP30_filter
    awk '{if ($6>=0.2 && $5>=25 && $7<=30 && $7>=5) {print}}' "$id" | wc -l \
        >> count/maf0.2_Q25_DP30_filter
done

# ---------------- 9. Count variant types -----------------------
echo "Counting variant types..."
ls *_specific | while read id; do
    base=${id%%.*}
    for j in SNP MNP DEL INS; do
        echo -ne "${base}\t${j}\t" >> count/type_count_0.2_Q25_DP30_homo
        awk -v a="${j}" '{if ($6>=0.2 && $5>=25 && $7<=30 && $7>=5 && $10==a) {print}}' "$id" | wc -l \
            >> count/type_count_0.2_Q25_DP30_homo
    done
done

# ---------------- 10. SNP subclassification --------------------
echo "Counting six SNP substitution types..."
for id in *_specific; do
    base=${id%%.*}
    awk '{if ($10=="SNP" && $6>=0.2 && $5>=25 && $7<=30 && $7>=5) {print $0}}' "$id" \
        > SNP_type/${base}.snp_0.2_Q25_DP30_homo
done

for id in *.snp_0.2_Q25_DP30_homo; do
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
    }' "$id" >> SNP_count/SNP_change_count_0.2_Q25_DP30_homo
done

# ---------------- 11. Annotation matching ---------------------
echo "Annotating variants..."
ls *_specific | while read id; do
    base=${id%%.*}
    awk -v FS='\t' '{if ($10=="DEL" || $10=="INS"){print $0"\t"$2+1"\t"$10} else if ($10=="SNP"){print $0"\t"$2"\t""SNV"} else{print $0"\t"$2"\t""MNV"}}' "$id" \
        > anno/${base}.corr
done

ls *.corr | while read id; do
    base=${id%%.*}
    awk -v FS='\t' 'FNR==NR {a[$3"_"$4"_"$6]=$7"\t"$8"\t"$9"\t"$10; next} ($1"_"$11"_"$12 in a){print $0"\t"a[$1"_"$11"_"$12]}' \
        "$ANNOTATION_DIR/${base}.decompose_normalize_annot.tsv" "$id" \
        > ${base}.anno
done

# ---------------- 12. Count annotation categories --------------
echo "Counting annotation summaries..."
ls *.anno | while read id; do
    base=${id%%.*}
    perl anno_genebody_AF02_Q25_DP30.pl "$id" anno_count/${base}.3.AF0.2_Q25_DP30.count_homo
    awk -v id=${base} '{print id "\t" $0}' anno_count/${base}.3.AF0.2_Q25_DP30.count_homo \
        >> anno_count/count/3_count_AF0.2_Q25_DP30_homo
done

ls *.anno | while read id; do
    base=${id%%.*}
    perl anno_CDS.pl "$id" anno_count/${base}.count_AF0.2_Q25_DP30_homo
    awk -v id=${base} '{print id "\t" $0}' anno_count/${base}.count_AF0.2_Q25_DP30_homo \
        >> anno_count/count/CDS_count_AF0.2_Q25_DP30_homo
done


