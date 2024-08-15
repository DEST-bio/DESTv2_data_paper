	gowinda_dir=$HOME/Dropbox/stockholm/projects/droseu/DEST/gowinda
	d_mel_ref_dir=$HOME/Dropbox/stockholm/projects/d_mel_ref
	snps_dir=$HOME/Dropbox/stockholm/projects/droseu/DEST/rdata/dest2_seasonality
	dest_dir=$HOME/Dropbox/stockholm/projects/droseu/DEST/vcf


Outlier window regions come from processing in R: 


Obtain the SNPs within the outlier windows.

	intersectBed -a ./dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz -b xtx_outlier_windows.bed > ./xtx_outlier_windows_snps.vcf

Collect just the chromosome and position information for all SNPs

	cut -f 1-2 ./xtx_outlier_windows_snps.vcf > ./xtx_outlier_windows_candidate_snps.tab

Get a file for all SNPs

	zcat ./dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz | grep -v "^#" | cut -f 1-2 > ./all_snps.tab

Use this snps table as input to GOwinda

XTX

	java -jar ${gowinda_dir}/Gowinda-1.12.jar \
		--snp-file ${dest_dir}/all_snps.tab \
		--candidate-snp-file ${snps_dir}/xtx_outlier_windows_candidate_snps.tab \
		--gene-set-file ${gowinda_dir}/funcassociate_go_associations.txt \
		--annotation-file ${d_mel_ref_dir}/dmel-all-r6.12.gtf \
		--simulations  100000 \
		--min-significance 1 \
		--gene-definition gene \
		--threads 4 \
		--output-file ${snps_dir}/xtx_outliers_gowinda_goea_results.txt \
		--mode gene \
		--min-genes 5


C2

	java -jar ${gowinda_dir}/Gowinda-1.12.jar \
		--snp-file ${dest_dir}/all_snps.tab \
		--candidate-snp-file ${snps_dir}/c2_outlier_windows_candidate_snps.tab \
		--gene-set-file ${gowinda_dir}/funcassociate_go_associations.txt \
		--annotation-file ${d_mel_ref_dir}/dmel-all-r6.12.gtf \
		--simulations  100000 \
		--min-significance 1 \
		--gene-definition gene \
		--threads 4 \
		--output-file ${snps_dir}/c2_outliers_gowinda_goea_results.txt \
		--mode gene \
		--min-genes 5


Extract only terms with FDR < 0.05 from results

	awk -F'\t' ' $5 < 0.05 ' ${snps_dir}/c2_outliers_gowinda_goea_results.txt > ${snps_dir}/c2_outliers_gowinda_goea_results.fdr0.05.txt

Extract the set of genes found near candidate SNPs by GOwinda
	awk 'BEGIN{OFS=FS="\t"}{print $10}' xtx_outliers_gowinda_goea_results.txt | sed 's;,;\n;g' | sort -k1,1 | uniq > xtx_outlier_gowinda_genes.txt
	awk 'BEGIN{OFS=FS="\t"}{print $10}' c2_outliers_gowinda_goea_results.txt | sed 's;,;\n;g' | sort -k1,1 | uniq > c2_outlier_gowinda_genes.txt


Machado et al. 2021

	sed -i '/^#/d' dmel-all-r5.5.gff

	grep -P "\tgene\t" dmel-all-r5.5.gff > dmel-all-genes-r5.5.gff

	intersectBed -b machado_2021_seasonalsnps.bed -a dmel-all-genes-r5.5.gff | awk 'BEGIN{FS=OFS="\t"}{print $9}' | sed 's/ID=\(FBgn[[:digit:]]\+\).*$/\1/g' | sort -k1,1 | uniq > machado_2021_seasonalsnps_genes.txt

Nunez et al. 2024

	d_mel_ref_dir=$HOME/Dropbox/stockholm/projects/d_mel_ref

	intersectBed -b outlier_windows.bed -a ${d_mel_ref_dir}/dmel-all-genes-r6.12.gff | awk 'BEGIN{FS=OFS="\t"}{print $9}' | sed 's/ID=\(FBgn[[:digit:]]\+\).*$/\1/g' | sort -k1,1 | uniq > nunes_2024_outlier_windows_genes.txt

DEST v2 Fgt

	awk 'BEGIN{FS=OFS="\t"}{if ( $3 > 0.2 ) print $4,$5,$6}' hierFST_k4_10.sort.csv > hierFST_k4_10.sort_fgt0.2.bed

	intersectBed -a ${dest_dir}/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz -b hierFST_k4_10.sort_fgt0.2.bed > ./hierFST_k4_10.sort_fgt0.2_snps.vcf

Collect just the chromosome and position information for all SNPs

	cut -f 1-2 ./hierFST_k4_10.sort_fgt0.2_snps.vcf > ./hierFST_k4_10.sort_fgt0.2_snps.tab

	intersectBed -b hierFST_k4_10.sort_fgt0.2_snps.tab -a ${d_mel_ref_dir}/dmel-all-genes-r6.12.gff | awk 'BEGIN{FS=OFS="\t"}{print $9}' | sed 's/ID=\(FBgn[[:digit:]]\+\).*$/\1/g' | sort -k1,1 | uniq > hierFST_k4_10.sort_fgt0.2_genes.txt