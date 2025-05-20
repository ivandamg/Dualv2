# NCR_Nathalie



## 1. Download data

Create a file containing address of files to download

          mkdir 01_RawReads
          echo "" > 01_RawReads/list_to_download # create empty file and copy paste links to data

Download with wget

          sbatch --partition=pshort_el8 --job-name=wget --time=0-02:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=wget.log --error=wget.err --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/30_DualRNaseq2/01_RawData; wget -i list_to_download"

## 2. Change filename to discard unique names

      cp 127_N_L1_R1_001_WVXjl0uZX5GP.fastq.gz NCR127_R1.fastq.gz
      cp 127_N_L1_R2_001_YcsxCGcFiFBs.fastq.gz NCR127_R2.fastq.gz
      cp 16_N_L1_R1_001_bF6zPIAMRHym.fastq.gz NCR16_R1.fastq.gz
      cp 16_N_L1_R2_001_Dpr0Mk78eRzm.fastq.gz NCR16_R2.fastq.gz
      .
      .
      .
      mkdir ../02_NewName
      mv N* ../02_NewName

## 3. Extract UMI

    mkdir ../03_UMIExtract
    for FILE in $(ls 12*_R1.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)_UMI --time=0-08:00:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=UMI_$(echo $FILE | cut -d'_' -f1).out --error=UMI_$(echo $FILE | cut -d'_' -f1).error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/02_NewName ;module load UMI-tools/1.0.1-foss-2021a; umi_tools extract --bc-pattern=NNNNNNNNNNNN --stdin=$FILE --stdout=../03_UMIExtract/$(echo $FILE | cut -d'_' -f1)_extracted_R1.fastq.gz --read2-in=$(echo $FILE | cut -d'_' -f1)_R2.fastq.gz --read2-out=../03_UMIExtract/$(echo $FILE | cut -d'_' -f1)_extracted_R2.fastq.gz"; sleep  1; done

## 4. Check quality

    for FILE in $(ls *.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastQC --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=$(echo $FILE | cut -d'_' -f1)_fastQC.out --error=$(echo $FILE | cut -d'_' -f1)_fastQC.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/03_UMIExtract ; module load FastQC; fastqc -t 4 $FILE"; sleep  1; done


## 5. fastp trimming and filter reads

        for FILE in $(ls *extracted_R1.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/03_UMIExtract ; module load FastQC; module load fastp/0.23.4-GCC-10.3.0; fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_extracted_R2.fastq.gz --out1 ../04_TrimmedData2/$(echo $FILE | cut -d'_' -f1)_R1_trimmed.fastq.gz --out2 ../04_TrimmedData2/$(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz -h ../04_TrimmedData2/$(echo $FILE | cut -d',' -f1)_fastp.html --detect_adapter_for_pe --trim_poly_g --length_required 99 --thread 4; fastqc -t 4 ../04_TrimmedData2/$(echo $FILE | cut -d'_' -f1)_R1_trimmed.fastq.gz; fastqc -t 4 ../04_TrimmedData2/$(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz"; sleep  1; done


# 6. index genome Rhizobia

        sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/00_References; module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/00_References --genomeFastaFiles FribourgSMeliloti_Prokka.fna --sjdbGTFfile FribourgSMeliloti_Prokka_v2.gff --sjdbGTFfeatureExon CDS --sjdbOverhang 99 --genomeSAindexNbases 10"




# 7. Map reads to rhizobia

           mkdir ../05_MappedRhizobia/
           for FILE in $(ls *_R1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)_1STAR --time=0-12:00:00 --mem-per-cpu=16G --ntasks=4 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_STAR.out --error=$(echo $FILE | cut -d'_' -f1)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/04_TrimmedData2; STAR --runThreadN 8 --genomeDir /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/00_References --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix $(echo $FILE | cut -d'_' -f1)_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5919206202 --limitOutSJcollapsed 5000000; mv $(echo $FILE | cut -d'_' -f1)*.bam ../05_MappedRhizobia/"; sleep  1; done


# 8. Deduplicate

                    mkdir ../06_Deduplicated/
                        for FILE in $(ls *.bam); do echo $FILE; sbatch --partition=pibu_el8 --job-name=dedup_$(echo $FILE | cut -d'_' -f1) --time=0-22:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_UMIdedup.out --error=$(echo $FILE | cut -d'_' -f1)_UMIdedup.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/05_MappedRhizobia; module load UMI-tools/1.0.1-foss-2021a; module load SAMtools/1.13-GCC-10.3.0; samtools index $FILE; umi_tools dedup --extract-umi-method=read_id --stdin=$FILE --stdout=$(echo $FILE | cut -d'_' -f1)_dedup.bam --paired --unpaired-reads=discard --chimeric-pairs=discard --log=$(echo $FILE | cut -d'_' -f1)_dedup.log; mv $(echo $FILE | cut -d'_' -f1)_dedup.* ../06_Deduplicated/"; done
#samtools sort -o $(echo $FILE | cut -d'.' -f1)_sorted.bam $FILE ; samtools index $(echo $FILE | cut -d'_' -f1)_MappedMedicago_sorted.bam;

# 8. Counts

        for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1).out --error=FC_$(echo $FILE | cut -d'_' -f1).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p --countReadPairs -t CDS -g ID -a /data/projects/p495_SinorhizobiumMeliloti/20_Nathalie/00_References/FribourgSMeliloti_Prokka_v2.gff  -o CountsTableRhizobia_$(echo $FILE | cut -d'_' -f1).txt $FILE -T 8"; sleep  1; done





# 8. divide gff into proteins and others


        cat GCF_037023865.1_ASM3702386v1_genomic.gff |grep "gene_biotype=p" > Rhizobium_proteins.gff

        cat GCF_037023865.1_ASM3702386v1_genomic.gff |grep -v "gene_biotype=p" > Rhizobium_others.gff



# 10. Counts also multiple mapped reads once
###  10a. count  reads to others tRNA and rRNA

        for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia/Rhizobium_others.gff  -o CountsTableRhizobia_UniqueMultiple_Others_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done

###  10b.  count  reads to proteins

       for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia/Rhizobium_proteins.gff  -o CountsTableRhizobia_UniqueMultiple_Proteins_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done

# 11. IGView viz create bai file for vizualisation

       for FILE in $(ls *.bam); do echo $FILE; sbatch --partition=pshort_el8 --job-name=index_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_index.out --error=$(echo $FILE | cut -d'_' -f1,2)_index.error --mail-type=END,FAIL --wrap "module load SAMtools/1.13-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/03_TrimmedData; samtools index $FILE"; done 

