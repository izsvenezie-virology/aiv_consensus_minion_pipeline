#USAGE: bash minion_aiv_consensus.sh READS OUTPUT_DIR NUMBER_CPU NAME
#READS: single fastq.gz file with minion data for which create consensus sequence
#OUTPUT_DIR: output directory in which do the analysis and put results (should be a novel one to avoid any issue)
#NUMBER_CPU: maximum number of cpu used by this script
#NAME: optional; add a string here to have it appear in fasta heander consensus file

#VARIABLES
export data=${1}
export out=${2}
export N_cpu=${3}
export name=${4}
export source_dir=$(dirname "$0")
export db=${source_dir}/data/AIV_complete_sequences_reduced.fa

#DATA FILTERING
if [ ! -d ${out} ]; then mkdir ${out}; fi;
zcat ${data} | chopper --minlength 200 --maxlength 2500 --quality 10 2>${out}/file.fil.err | gzip >${out}/file.fil.fastq.gz

#FIND REFERENCE
zcat ${out}/file.fil.fastq.gz | awk '{getline; print $0; getline; getline}' | fold -w 100 | awk '{if (length($0) >= 80) {print}}' | awk 'BEGIN {i=1} {print ">seq"i"__"length($0)"\n"$0; i++}' >${out}/file.fil.fa
if [ ! -e ${db}.nsq ]; then gzip -d ${db}.gz; makeblastdb -dbtype nucl -in ${db} >${db}.log 2>${db}.err; fi
blastn -db ${db} -query ${out}/file.fil.fa -outfmt '6 std stitle' -evalue 0.0000000001 -num_threads ${N_cpu} 2>${out}/file.tab.log | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$13}' | sed 's/^.*__\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\$/\1\t\2\t\3\t\4\t\5/' | awk -F "\t" '{if ( ($4>=$1*0.9) && ($3>=90) ) {print $2"\t"$5}}' | sed 's/ /_/g' | perl ${source_dir}/src/filtra_blast_tab.pl ${out}/file_mod.fil.tab.entry_title2hits ${out}/file_mod.fil.tab.entry_hit ${out}/file_mod.fil.tab.entry_title >${out}/file_mod.fil.log 2>${out}/file_mod.fil.err
bash ${source_dir}/src/estrai_best_hit_influenza_v2.0.sh ${out}/file_mod.fil.tab.entry_title 0 | perl ${source_dir}/src/accnum2sequence.pl ${out}/file_mod.fil.tab ${db} | perl ${source_dir}/src/remove_degeneration.pl >${out}/reference.fa 2>${out}/reference.fa.err

#READ ALIGNMENT
minimap2 -a -t ${N_cpu} -x map-ont ${out}/reference.fa ${out}/file.fil.fastq.gz 2>${out}/alignment_minimap2.err | samtools view -bh - 2>${out}/alignment_sam2bam.err | samtools sort -T ${out}/alignment_nosorting_temp --threads ${N_cpu} -o ${out}/alignment_possorting.bam - >${out}/alignment_possorting.bam.log 2>${out}/alignment_possorting.bam.err

#CONSENSUS CREATION
samtools view -F 3844 ${out}/alignment_possorting.bam 2>/dev/null | cut -f 3,4,6,10 | perl ${source_dir}/src/controllo_sito_clivaggio_AIV.pl ${out}/reference.fa >${out}/indel_cleavage_site.txt 2>${out}/indel_cleavage_site.err
samtools view -F 3844 ${out}/alignment_possorting.bam 2>/dev/null | cut -f 3,4,6,10,11 | perl ${source_dir}/src/create_consensus.pl ${out}/reference.fa 10 100 10 ${out}/degeneration.txt ${name} 2>${out}/consensus_raw.fa.log | perl ${source_dir}/src/aggiunta_sito_clivaggio_AIV.pl ${out}/indel_cleavage_site.txt >${out}/consensus.fa 2>${out}/consensus.fa.log
