# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:42:11 2021

ver 0.31 - "Enjoy exaggerate" (Public beta) 11.21.2021
    _x_ improved rarefaction plots

Purpose: 
    1. generate a bash script for the preprocessing and alignment of 
RNAseq data with UMIs.
    2. QC RNAseq data for duplicates using UMI and rRNA contamination
    3. Additional QC using UMI aware sequence rarefaction
    4. Calculate read coverage for a given feature file
    5. Combine coverage into a single table for downstream (ie DESeq2) analysis

Control File format:
    Every control file must have the following information in the following format
    
    #data_type	variable	value
    meta	set_name	<string>
    meta	data_dir	<path to fastq files>
    meta	work_dir	<path to a 'temp' directory>
    meta	output_dir	<path to a results directory>
    meta	genome_fa	<path to a reference fasta file>
    meta	genome_features	<path to a gff file>
    meta	intron_max	<integer>
    meta	adapter_seq_R1	<string>
    meta	adapter_seq_R2	<string>
    meta	prefix	<string>
    meta	fastq_name_template	<string>
    meta	qc_locus	<sting>
    sample	<string>
    sample	<string
    
Command format:
    # generate Align script 
    python windchime.py -a -i ctrl_file.tab -o run_windchime_star.sh
    # run rarefaction
    python windchime.py -r -i ctrl_file.tab
    # evaluate qc
    python windchime.py -e -i ctrl_file.tab -o evaluate_run.txt
    # generate coverage script
    python windchime.py -c -i ctrl_file.tab -o coverage_run.sh
    # combine coverage data
    python windchime.py -t -i ctrl_file.tab -o coverage_table.txt

#

@author: Pieter Spealman
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

# io
parser.add_argument('-i', '--inputfile_name')
parser.add_argument('-o', '--outfile_name')

# options
parser.add_argument('-mem', '--mem_in_GB')
parser.add_argument('-tol', '--set_tolerance')
parser.add_argument('-rep', '--use_replicates', action='store_true')


# functions
parser.add_argument('-r', '--rarefaction', action='store_true')
parser.add_argument('-a', '--align', action='store_true')
parser.add_argument('-e', '--evaluate', action='store_true')
parser.add_argument('-c', '--coverage', action='store_true')
parser.add_argument('-t', '--table_coverage_build', action='store_true')

args = parser.parse_args()  

# Handle default parameters
if not args.mem_in_GB:
    mem = 60
else:
    mem = int(args.mem_in_GB)
    
if not args.set_tolerance:
    tol = 0.500
else:
    tol = float(args.set_tolerance)
    if tol > 1:
        tol /= 100
        
#

def load_cmd_file():
    if args.inputfile_name:
        command_file_name = args.inputfile_name
        command_file = open(command_file_name)
        
        if args.use_replicates:
            cmd_dict = {
                'set_name':'',
                'data_dir':'',
                'work_dir':'',
                'output_dir':'',
                'genome_fa':'',
                'genome_features':'',
                'intron_max':'',
                'adapter_seq_R1':'',
                'adapter_seq_R2':'',
                'fastq_name_template':'',
                'prefix':'',
                'qc_locus':'',
                'samples':{},
                }
            
            for line in command_file:
                if line[0] != '#':
                    line = line.strip()
                
                    if line.split('\t')[0] == 'meta':
                        meta_val = line.split('\t')[1]
                        cmd_dict[meta_val] = line.split('\t')[2]
                        
                    if line.split('\t')[0] == 'sample':
                        strain = line.split('\t')[1]
                        replicate = line.split('\t')[2]
                        
                        if strain not in cmd_dict['samples']:
                            cmd_dict['samples'][strain] = set()
                            
                        cmd_dict['samples'][strain].add(replicate)
                        
        if not args.use_replicates:
            cmd_dict = {
                'set_name':'',
                'data_dir':'',
                'work_dir':'',
                'output_dir':'',
                'genome_fa':'',
                'genome_features':'',
                'intron_max':'',
                'adapter_seq_R1':'',
                'adapter_seq_R2':'',
                'fastq_name_template':'',
                'prefix':'',
                'qc_locus':'',
                'samples':set(),
                }
            
            for line in command_file:
                if line[0] != '#':
                    line = line.strip()
                
                    if line.split('\t')[0] == 'meta':
                        meta_val = line.split('\t')[1]
                        cmd_dict[meta_val] = line.split('\t')[2]
                        
                    if line.split('\t')[0] == 'sample':
                        strain = line.split('\t')[1]
                                                
                        cmd_dict['samples'].add(strain)
                                        
        command_file.close()
            
        process = True
        for each in cmd_dict:
            if cmd_dict[each] == '':
                outline = ('Please enter a value for {}').format(each)
                print(outline)
                process = False
        
        return(cmd_dict, process)
    
    else:
        print('Please specify a input control file')
            
def make_text(cmd_dict):
    outfile_name = args.outfile_name
    outfile = open(outfile_name, 'w')
    
    set_name = cmd_dict['set_name']
    qc_locus = cmd_dict['qc_locus']
    qc_text = qc_locus.replace(':','_')
        
    outline = ('#!/bin/bash\n'
               '#\n'
               '#SBATCH --verbose\n'
               '#SBATCH --job-name={set_name_py}_windchime\n'
               '#SBATCH --output={set_name_py}_windchime_%j.out\n'
               '#SBATCH --error={set_name_py}_windchime_%j.err\n'
               '#SBATCH --time=96:00:00\n'
               '#SBATCH --nodes=4\n'
               '#SBATCH --mem={mem}GB\n'
               '# 11.04.21 Pieter Spealman \n'
               '#=============================\n'
               '# 0. Load Modules\n'
               '#=============================\n'
               '\tmodule load star/intel/2.7.6a\n'
               '\tmodule load samtools/intel/1.12\n'
               '\tmodule load cutadapt/3.1\n'
               '#=============================\n'
               '# 1. Set Variables, Directories\n'
               '#=============================\n'
               '# Set name will form the part of the namespace and directory structure. \n'
               '\tset_name={set_name_py}\n'
               '\tadapter_seq_R1={adapter_seq_R1_py}\n'
               '\tadapter_seq_R2={adapter_seq_R2_py}\n'
               '\tprefix={prefix_py}\n'
               '\tdata_dir={data_dir_py}\n'
               '\techo "Processing" $set_name from $data_dir "..."\n'
               '\twork_dir={work_dir_py}/STAR_$set_name/\n'
               '\tfastq_dir=$work_dir/fastq/\n'
               '\tstar_idx_dir=$work_dir/idx/\n'
               '\trrna_idx=$star_idx_dirrrna/\n'
               '\tgenome_idx=$star_idx_dir/genome/\n'
               '\ttmp_dir=$work_dir/output/\n'
               '\tQC_dir=$tmp_dirQC/\n'
               '\tProcessed_dir={output_dir_py}\n'
               '#=============================\n'
               '# 2. Make directories\n'
               '#=============================\n'
               '\tmkdir -p $work_dir\n'
               '\tmkdir -p $data_dir\n'
               '\tmkdir -p $fastq_dir\n'
               '\tmkdir -p $star_idx_dir\n'
               '\tmkdir -p $rrna_idx\n'
               '\tmkdir -p $genome_idx\n'
               '\tmkdir -p $tmp_dir\n'
               '\tmkdir -p $QC_dir\n'
               '\tmkdir -p $Processed_dir\n'
               '#=============================\n'
               '# 3. Set STAR variables\n'
               '#=============================\n'
               '\tcp {genome_fa_py} $genome_idx/reference_genome.fa\n'
               '\tgenome_fa=$genome_idx/reference_genome.fa\n'
               '\tnproc=8\n'
               '\tnmismatch=2\n'
               '\tQC_suffix=_genome_Aligned.out\n'
               '\t# set STAR alignments\n'
               '\tnc_align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 100 --outFilterMismatchNmax $nmismatch  --alignIntronMax {intron_max}"\n'
               '\talign_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 1 --outFilterScoreMin 20 --outFilterMultimapNmax 3 --outFilterMatchNmin 20 --outFilterMismatchNmax $nmismatch --outFilterIntronMotifs RemoveNoncanonical --scoreGap -8 --scoreGapNoncan -16 --alignIntronMax {intron_max}"\n'
               '\t# set output format\n'
               '\tSAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMprimaryFlag AllBestScore"\n'
               '\t#\n'
               '\techo "building genome index..."\n'
               '\tSTAR --runThreadN $nproc --runMode genomeGenerate --genomeDir $genome_idx --genomeFastaFiles {genome_fa_py} --genomeSAindexNbases 4\n'
               '#=============================\n'
               '# 4. Preprocess samples\n'
               '#=============================\n').format(
                   set_name_py=set_name,
                   mem = mem,
                   adapter_seq_R1_py=cmd_dict['adapter_seq_R1'],
                   adapter_seq_R2_py=cmd_dict['adapter_seq_R2'],
                   prefix_py=cmd_dict['prefix'],
                   data_dir_py=cmd_dict['data_dir'],
                   work_dir_py=cmd_dict['work_dir'],
                   output_dir_py=cmd_dict['output_dir'],
                   genome_fa_py=cmd_dict['genome_fa'],
                   intron_max = int(cmd_dict['intron_max'])
                   )
                   
    outfile.write(outline)
    
    icpy = 0
    
    icpy_assigned_dict = {} 
    
    for strain_py in cmd_dict['samples']:
        if args.use_replicates:
            for replicate_py in cmd_dict['samples'][strain_py]:
                if icpy not in icpy_assigned_dict:
                    icpy_assigned_dict[icpy] = {'strain':strain_py,'replicate':replicate_py}
                else:
                    print('Error: icpy duplicate generated')
                    
                input_file_R1_is = cmd_dict['fastq_name_template']
                input_file_R2_is = cmd_dict['fastq_name_template'].replace('n01', 'n02')
                    
                outline = ('#Preprocess_{icpy}\n'
                       '\tstrain={strain_py}\n'
                       '\treplicate={replicate_py}\n'
                       '\t#\n'
                       '\t\tsample=$set_name-$strain-$replicate\n'
                       '\t\tinput_file_R1={input_file_R1_py}\n'
                       '\t\toutput_file_R1={output_file_R1_py}\n'
                       '\t\tinput_file_R2={input_file_R2_py}\n'
                       '\t\toutput_file_R2={output_file_R2_py}\n'
                       '\t#\n'
                       '\t\toutputfile_{icpy}=$sample\n'
                       '\t\techo "Starting on " $sample\n'
                       'echo $input_file_R1.fastq.gz\n'
                       '\t\tcp $data_dir$input_file_R1.fastq.gz $fastq_dir$output_file_R1.fastq.gz\n'
                       '\t\tgunzip -f $fastq_dir$output_file_R1.fastq.gz\n'
                       '\t\tcp $data_dir$input_file_R2.fastq.gz $fastq_dir$output_file_R2.fastq.gz\n'
                       '\t\tgunzip -f $fastq_dir$output_file_R2.fastq.gz\n'
                       '\t# cutadapt #\n'
                       '\t\tcutadapt --cores=20 -e 0.12 -m 30 -a $adapter_seq_R1 -A $adapter_seq_R2 -o $fastq_dir/{trimmed_R1} -p $fastq_dir/{trimmed_R2} $fastq_dir/$output_file_R1.fastq $fastq_dir/$output_file_R2.fastq\n'
                       '\t# umi extract #\n'
                       '\t\tumi_tools extract --either-read-resolve=quality --bc-pattern=NNNNNNNNNNNN --bc-pattern2=NNNNNNNNNNNN -I $fastq_dir/{trimmed_R1} --read2-in=$fastq_dir/{trimmed_R2} --stdout=$fastq_dir/{umi_R1} --read2-out=$fastq_dir/{umi_R2} --log=$fastq_dir/$sample_processed.log\n'
                       '\t# set sample star variables #\n'
                       '\t\triboseq_fa_{icpy}_R1=$fastq_dir/{umi_R1}\n'
                       '\t\triboseq_fa_{icpy}_R2=$fastq_dir/{umi_R2}\n'
                       '\t\tgriboprefix_{icpy}=$tmp_dir/{sample_genome}\n'
                       '\t# Clean up #\n'
                       '\t\trm $fastq_dir$output_file_R1.fastq\n'
                       '\t\trm $fastq_dir$output_file_R2.fastq\n'
                       '\t\trm $fastq_dir/{trimmed_R1}\n'
                       '\t\trm $fastq_dir/{trimmed_R2}\n').format(
                           icpy=icpy,
                           strain_py=strain_py,
                           replicate_py=replicate_py,
                           input_file_R1_py=input_file_R1_is,
                           input_file_R2_py=input_file_R2_is,
                           output_file_R1_py="${sample}_1",
                           output_file_R2_py="${sample}_2",
                           trimmed_R1="${output_file_R1}_trimmed.fastq",
                           trimmed_R2="${output_file_R2}_trimmed.fastq",
                           umi_R1="${output_file_R1}_umi.fastq",
                           umi_R2="${output_file_R2}_umi.fastq",
                           sample_genome="${sample}_genome_"
                           )
                           
                outfile.write(outline)
            
                icpy+=1
                           
        if not args.use_replicates:
            if icpy not in icpy_assigned_dict:
                icpy_assigned_dict[icpy] = {'strain':strain_py}
            else:
                print('Error: icpy duplicate generated')
                
            if not args.use_replicates:
                input_file_R1_is = cmd_dict['fastq_name_template']
                input_file_R2_is = cmd_dict['fastq_name_template'].replace('n01', 'n02')
                
            outline = ('#Preprocess_{icpy}\n'
                       '\tstrain={strain_py}\n'
                       '\t\n'
                       '\t#\n'
                       '\t\tsample=$set_name-$strain\n'
                       '\t\tinput_file_R1={input_file_R1_py}\n'
                       '\t\toutput_file_R1={output_file_R1_py}\n'
                       '\t\tinput_file_R2={input_file_R2_py}\n'
                       '\t\toutput_file_R2={output_file_R2_py}\n'
                       '\t#\n'
                       '\t\toutputfile_{icpy}=$sample\n'
                       '\t\techo "Starting on " $sample\n'
                       '\t\techo $input_file_R1.fastq.gz\n'
                       '\t\tcp $data_dir$input_file_R1.fastq.gz $fastq_dir$output_file_R1.fastq.gz\n'
                       '\t\tgunzip -f $fastq_dir$output_file_R1.fastq.gz\n'
                       '\t\tcp $data_dir$input_file_R2.fastq.gz $fastq_dir$output_file_R2.fastq.gz\n'
                       '\t\tgunzip -f $fastq_dir$output_file_R2.fastq.gz\n'
                       '\t# cutadapt #\n'
                       '\t\tcutadapt --cores=20 -e 0.12 -m 30 -a $adapter_seq_R1 -A $adapter_seq_R2 -o $fastq_dir/{trimmed_R1} -p $fastq_dir/{trimmed_R2} $fastq_dir/$output_file_R1.fastq $fastq_dir/$output_file_R2.fastq\n'
                       '\t# umi extract #\n'
                       '\t\tumi_tools extract --either-read-resolve=quality --bc-pattern=NNNNNNNNNNNN --bc-pattern2=NNNNNNNNNNNN -I $fastq_dir/{trimmed_R1} --read2-in=$fastq_dir/{trimmed_R2} --stdout=$fastq_dir/{umi_R1} --read2-out=$fastq_dir/{umi_R2} --log=$fastq_dir/$sample_processed.log\n'
                       '\t# set sample star variables #\n'
                       '\t\triboseq_fa_{icpy}_R1=$fastq_dir/{umi_R1}\n'
                       '\t\triboseq_fa_{icpy}_R2=$fastq_dir/{umi_R2}\n'
                       '\t\tgriboprefix_{icpy}=$tmp_dir/{sample_genome}\n'
                       '\t# Clean up #\n'
                       '\t\trm $fastq_dir$output_file_R1.fastq\n'
                       '\t\trm $fastq_dir$output_file_R2.fastq\n'
                       '\t\trm $fastq_dir/{trimmed_R1}\n'
                       '\t\trm $fastq_dir/{trimmed_R2}\n').format(
                           icpy=icpy,
                           strain_py=strain_py,
                           input_file_R1_py=input_file_R1_is,
                           input_file_R2_py=input_file_R2_is,
                           output_file_R1_py="${sample}_1",
                           output_file_R2_py="${sample}_2",
                           trimmed_R1="${output_file_R1}_trimmed.fastq",
                           trimmed_R2="${output_file_R2}_trimmed.fastq",
                           umi_R1="${output_file_R1}_umi.fastq",
                           umi_R2="${output_file_R2}_umi.fastq",
                           sample_genome="${sample}_genome_"
                           )

            outfile.write(outline)
            
            icpy+=1
            
    outline = ('#=============================\n'
               '# 5. Run STAR\n'
               '#=============================\n')
    outfile.write(outline)
            
    for icpy in icpy_assigned_dict:
        outline = ('\tSTAR --runThreadN $nproc '
                   '--genomeDir $genome_idx '
                   '--readFilesIn $riboseq_fa_{icpy}_R1 $riboseq_fa_{icpy}_R2 '
                   '--outFileNamePrefix $griboprefix_{icpy} '
                   '$SAM_params $align_params\n').format(icpy=icpy)
        
        outfile.write(outline)
        
    outline = ('#=============================\n'
           '# 6. Run QC\n'
           '#=============================\n'
           '\techo "Beginning" $set_name "bam sorting and indexing ..."\n')
    outfile.write(outline)
    
    for icpy in icpy_assigned_dict:
        outline = ('\t\tQC_name=$outputfile_{icpy}\n'
                   '\t\t# prepare bams #\n'
                   '\t\tsamtools sort -m 5000000000 {temp_QCname_py}$QC_suffix.bam -o {temp_QCname_py}.sorted.bam\n'
                   '\t\tsamtools index {temp_QCname_py}.sorted.bam\n'
                   '\t#\n'
                   '\t\tumi_tools dedup -I {temp_QCname_py}.sorted.bam --output-stats={temp_QCname_py}.deduplicated -S {temp_QCname_py}.bam\n'
                   '\t#\n'
                   '\t\tsamtools index {temp_QCname_py}.bam\n'
                   '\t\tsamtools view -b -h {temp_QCname_py}.bam "{qc_locus}" > {temp_QCname_py}_{qc_text}_output.bam'
                   '\t# output QC #\n'
                   '\t\techo "total aligned:" > {temp_QCname_py}_total_stats.log\n'
                   '\t\tsamtools stats {temp_QCname_py}.sorted.bam >> {temp_QCname_py}_total_stats.log\n'
                   '\t\techo "total after deduplication of UMI:" > {temp_QCname_py}_dedup_stats.log\n'
                   '\t\tsamtools stats {temp_QCname_py}.bam >> {temp_QCname_py}_dedup_stats.log\n'
                   '\t\techo "total deduplication QC locus" > {temp_QCname_py}_{qc_text}_stats.log\n'
                   '\t\tsamtools stats {temp_QCname_py}_{qc_text}_output.bam >> {temp_QCname_py}_{qc_text}_stats.log\n'
                   '\t# Clean up #\n'
                   '\t\tmv {temp_QCname_py}.bam* $Processed_dir/\n'
                   '\t\tmv {temp_QCname_py}*.log $Processed_dir/\n'
                   '\t\tmv {temp_QCname_py}*.pdf $Processed_dir/\n').format(
                       icpy=icpy,
                       temp_QCname_py="${tmp_dir}/${QC_name}",
                       qc_locus=qc_locus,
                       qc_text=qc_text)
        outfile.write(outline)
    
    outfile.close()
   
def plot_rarefaction(x, y, ideal_y, halfline, outfile_name, sample_name):
        
    title_label = ('Rarefaction curve for {}').format(sample_name)
    
    axvh_label = ('{}% percent novel').format(tol*100)
    plt.figure()
    
    plt.plot(x, y, label='Observed rarefaction')
    plt.plot(x, ideal_y, 'r--', label='Ideal rarefaction')
    
    plt.title(title_label)
    plt.ylabel('Unique reads')
    plt.xlabel('Sequenced reads')
    
    hfy = halfline[1]
    hfx = halfline[0]

    plt.axhline(color='k', linewidth=0.5, y=hfy, label=axvh_label)
    plt.axvline(color='k', linewidth=0.5, x=hfx)
    
    plt.legend()
       
    plt.savefig(outfile_name + '.pdf')
    plt.close()
    plt.clf()
    
def run_rarefaction(cmd_dict):
    fastq_dir = ('{work_dir}/STAR_{set_name}/fastq').format(
        work_dir = cmd_dict['work_dir'],
        set_name = cmd_dict['set_name'])
    
    set_name = cmd_dict['set_name']
        
    if args.use_replicates:
        for strain_py in cmd_dict['samples']:
            for replicate_py in cmd_dict['samples'][strain_py]:
                outfile_name = ('{}/{}-{}').format(cmd_dict['output_dir'],strain_py, replicate_py)
                #outfile = open(outfile_name+'.txt', 'w')
                
                read_ct = 0
                new = 0
                
                collision_dict = {}
                
                x = []
                y = []
                ideal_y = []
                halfline = [0,0]
                
                print('Processing: ')
                
                sample_1 = ('{work_dir}/STAR_{set_name}/fastq/{set_name}-{strain}-{replicate}').format(
                    work_dir = cmd_dict['work_dir'],
                    set_name = cmd_dict['set_name'],
                    strain = cmd_dict['strain'],
                    replicate = cmd_dict['replicate'])
                
                print(sample_1)
                                
                infile_1 = open(sample_1)
                
                ct = 0
                
                for line in infile_1:
                    if (read_ct % 10000 == 0) and (ct == 0):
                        #outline = ('{}\t{}\n').format(read_ct, new)
                        #outfile.write(outline)
                        #
                        x.append(read_ct)
                        y.append(new)
                        ideal_y.append(read_ct)
                        
                        if halfline == [0, 0]:
                            if round(new/max(1,read_ct),3) == tol:
                                halfline = [read_ct, new]                       
                        
                    if line[0]=='@' and ct == 0:
                        read_ct += 1
                        umi = line.split(' ')[0].split('_')[1]
                        
                    if ct == 1:
                        seq = line.strip()
                        
                        if seq not in collision_dict:
                            collision_dict[seq] = set()
                            
                        if umi not in collision_dict[seq]:
                            new += 1
                            collision_dict[seq].add(umi)
                                                
                    ct+=1
                    
                    if ct >= 4:
                        umi=''
                        seq=''
                        ct = 0
                
                
                infile_1.close()
                #outfile.close()
                
                sample_name = ('{strain}-{replicate}').format(
                    strain = cmd_dict['strain'],
                    replicate = cmd_dict['replicate'])
                
                plot_rarefaction(x, y, ideal_y, halfline, outfile_name, sample_name)

                
    if not args.use_replicates:
        for strain_py in cmd_dict['samples']:
            outfile_name = ('{}/{}').format(cmd_dict['output_dir'],strain_py)
            #outfile = open(outfile_name+'.txt', 'w')
            
            read_ct = 0
            new = 0
            
            collision_dict = {}
            
            x = []
            y = []
            ideal_y = []
            halfline = [0,0]
            
            print('Processing: ')
            
            sample_1 = ('{fastq_dir}/{set_name}-{strain_py}_1_umi.fastq').format(
                set_name=set_name,
                fastq_dir=fastq_dir, 
                strain_py=strain_py)
            
            print(sample_1)
                        
            infile_1 = open(sample_1)
            
            ct = 0
            
            for line in infile_1:
                if (read_ct % 10000 == 0) and (ct == 0):
                    #outline = ('{}\t{}\n').format(read_ct, new)
                    #outfile.write(outline)
                    #
                    x.append(read_ct)
                    y.append(new)
                    ideal_y.append(read_ct)
                    
                    if halfline == [0, 0]:
                        if round(new/max(1,read_ct),3) == tol:
                            halfline = [read_ct, new]                           
                    
                if line[0]=='@' and ct == 0:
                    read_ct += 1
                    umi = line.split(' ')[0].split('_')[1]
                    
                if ct == 1:
                    seq = line.strip()
                                        
                    if seq not in collision_dict:
                        collision_dict[seq] = set()
                        
                    if umi not in collision_dict[seq]:
                        new += 1
                        collision_dict[seq].add(umi)
                                            
                ct+=1
                
                if ct >= 4:
                    umi=''
                    seq=''
                    ct = 0
            
            
            infile_1.close()
            #outfile.close()
            
            plot_rarefaction(x, y, ideal_y, halfline, outfile_name, strain_py)
                                         
def stats_parser(stat_file_name):
    print(stat_file_name)
    try:
        stats_file = open(stat_file_name)
    
        for line in stats_file:
            if 'reads mapped and paired:' in line:
                line = line.strip()
                
                return(int(line.split('\t')[2]))
    
    except:
        print("Can't open ", stat_file_name)
        return(False)
    
    
def summarize_stats(cmd_dict):
    
    outfile_name = args.outfile_name
    outfile = open(outfile_name, 'w')
        
    set_name = cmd_dict['set_name']
    output_dir = cmd_dict['output_dir']
    
    qc_locus = cmd_dict['qc_locus']
    qc_text = qc_locus.replace(':','_')
    
    if args.use_replicates:
        header = ('Strain\tReplicate'
                  '\tTotal_reads'
                  '\tUMI_unique_reads\tUMI_unique_pct'
                  '\tQC_locus_reads\tQC_locus_pct\n')
        print(header)
        outfile.write(header)
        
        for strain_py in cmd_dict['samples']:
            for replicate_py in cmd_dict['samples'][strain_py]:
                print(strain_py, replicate_py)
                
                total_file_name = ('{output_dir}/{set_name}-{strain_py}-{replicate_py}_total_stats.log').format(
                    output_dir=output_dir,
                    set_name = set_name,
                    strain_py = strain_py, 
                    replicate_py = replicate_py)
                rrna_file_name = ('{output_dir}/{set_name}-{strain_py}-{replicate_py}_{qc_text}_stats.log').format(
                    output_dir=output_dir,
                    set_name = set_name,
                    strain_py = strain_py, 
                    replicate_py = replicate_py,
                    qc_text=qc_text)
                dedup_file_name = ('{output_dir}/{set_name}-{strain_py}-{replicate_py}_dedup_stats.log').format(
                    output_dir=output_dir,
                    set_name = set_name,
                    strain_py = strain_py, 
                    replicate_py = replicate_py)
                
                total_reads_py = stats_parser(total_file_name)
                umi_reads_py = stats_parser(dedup_file_name)
                rrna_reads_py = stats_parser(rrna_file_name)
                
                
                outline = ('{strain}\t{replicate}'
                              '\t{total_reads}'
                              '\t{umi_reads}\t{umi_pct}'
                              '\t{rrna_reads}\t{rrna_pct}\n').format(strain=strain_py,
                                replicate = replicate_py,
                                total_reads = total_reads_py,
                                umi_reads = umi_reads_py, 
                                umi_pct = round(100*float(umi_reads_py/max(total_reads_py,1)),2),
                                rrna_reads = rrna_reads_py, 
                                rrna_pct = round(100*float(rrna_reads_py/max(umi_reads_py,1)),2),
                                )  
                print(outline)
                outfile.write(outline)
        
    if not args.use_replicates:
        header = ('Strain'
                  '\tTotal_reads'
                  '\tUMI_unique_reads\tUMI_unique_pct'
                  '\tQC_locus_reads\tQC_locus_pct\n')
        print(header)
        outfile.write(header)
        
        for strain_py in cmd_dict['samples']:
            print(strain_py)
            
            total_file_name = ('{output_dir}/{set_name}-{strain_py}_total_stats.log').format(
                output_dir=output_dir,
                set_name = set_name,
                strain_py = strain_py)
            rrna_file_name = ('{output_dir}/{set_name}-{strain_py}_{qc_text}_stats.log').format(
                output_dir=output_dir,
                set_name = set_name,
                strain_py = strain_py,
                qc_text=qc_text)
            dedup_file_name = ('{output_dir}/{set_name}-{strain_py}_dedup_stats.log').format(
                output_dir=output_dir,
                set_name = set_name,
                strain_py = strain_py)
            
            total_reads_py = stats_parser(total_file_name)
            umi_reads_py = stats_parser(dedup_file_name)
            rrna_reads_py = stats_parser(rrna_file_name)
            
            outline = ('{strain}'
                          '\t{total_reads}'
                          '\t{umi_reads}\t{umi_pct}'
                          '\t{rrna_reads}\t{rrna_pct}\n').format(strain=strain_py,
                            total_reads = total_reads_py,
                            umi_reads = umi_reads_py, 
                            umi_pct = round(100*float(umi_reads_py/max(total_reads_py,1)),2),
                            rrna_reads = rrna_reads_py, 
                            rrna_pct = round(100*float(rrna_reads_py/max(umi_reads_py,1)),2),
                            )  
            print(outline)
            outfile.write(outline)
            
    outfile.close()
    
def coverage_run(cmd_dict):
    
    outfile_name = args.outfile_name
    outfile = open(outfile_name, 'w')
    
    header = ('#!/bin/bash\n'
               '#\n'
               '#SBATCH --verbose\n'
               '#SBATCH --job-name={set_name_py}_windchime\n'
               '#SBATCH --output={set_name_py}_windchime_%j.out\n'
               '#SBATCH --error={set_name_py}_windchime_%j.err\n'
               '#SBATCH --time=96:00:00\n'
               '#SBATCH --nodes=1\n'
               '#SBATCH --mem={mem}GB\n'
               'module load bedtools/intel/2.29.2\n'
               'module load samtools/intel/1.12\n').format(
                   mem = mem,
                   set_name_py = cmd_dict['set_name'])
    outfile.write(header)
    
    genome_gtf = cmd_dict['genome_features']
    
    set_name = cmd_dict['set_name']
    
    output_dir_py=cmd_dict['output_dir']
    
    #qc_locus = cmd_dict['qc_locus']
    #qc_text = qc_locus.replace(':','_')
    
    if args.use_replicates:
        for strain_py in cmd_dict['samples']:
            for replicate_py in cmd_dict['samples'][strain_py]:
                print(strain_py, replicate_py)
                
                bam_file_name = ('{output_dir_py}/{set_name}-{strain_py}-{replicate_py}.bam').format(
                    output_dir_py = output_dir_py,
                    set_name = set_name,
                    strain_py = strain_py, 
                    replicate_py = replicate_py)
                
                tab_file_name = ('{output_dir_py}/{set_name}-{strain_py}-{replicate_py}.coverage.tab').format(
                    output_dir_py = output_dir_py,
                    set_name = set_name,
                    strain_py = strain_py, 
                    replicate_py = replicate_py)
    
                outline = ('bedtools coverage -counts -s -b {bam_file_name} '
            		'-a {genome_gtf} '
            		'> {tab_file_name}\n').format(
                    genome_gtf = genome_gtf,
                    bam_file_name = bam_file_name,
                    tab_file_name = tab_file_name)
                 
                #print(outline)
                outfile.write(outline)
                
    if not args.use_replicates:
        for strain_py in cmd_dict['samples']:
            print(strain_py)
            
            bam_file_name = ('{output_dir_py}/{set_name}-{strain_py}.bam').format(
                output_dir_py = output_dir_py,
                set_name = set_name,
                strain_py = strain_py)
            
            tab_file_name = ('{output_dir_py}/{set_name}-{strain_py}.coverage.tab').format(
                output_dir_py = output_dir_py,
                set_name = set_name,
                strain_py = strain_py)
                        
            outline = ('bedtools coverage -counts -s -b {bam_file_name} '
          		'-a {genome_gtf} '
          		'> {tab_file_name}\n').format(
                  genome_gtf = genome_gtf,
                  bam_file_name = bam_file_name,
                  tab_file_name = tab_file_name)
             
            #print(outline)
            outfile.write(outline)
                
    outfile.close()
    
def parse_coverage_file(tab_file_name, values_dict, gene_list, sample_uid, runmode):
    if runmode == 'gtf':
        tab_file = open(tab_file_name)
                
        for line in tab_file:
            if line[0] != '#':
                line = line.strip()
                
                #IV	sgd	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_name "COS7"; gene_source "sgd"; gene_biotype "protein_coding";	762
                if 'transcript' in line.split('\t')[2]:
                    gene_name = line.split('\t')[8]
                    gene_name = gene_name.split(';')[0]
                    gene_name = gene_name.split(' ')[1]
                    gene_name = gene_name.replace('"','')
                    
                    value = int(line.split('\t')[9])
                    
                    if gene_name not in values_dict:
                        values_dict[gene_name] = {}
                        gene_list.append(gene_name)
                        
                    if sample_uid not in values_dict[gene_name]:
                        values_dict[gene_name][sample_uid] = 0
                        
                    values_dict[gene_name][sample_uid] += value
        
        tab_file.close()
        
        return(values_dict, gene_list)
    
    if runmode == 'gff':
        tab_file = open(tab_file_name)
                
        for line in tab_file:
            if line[0] != '#':
                line = line.strip()
                
                #IV	sgd	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_name "COS7"; gene_source "sgd"; gene_biotype "protein_coding";	762
                if 'gene' in line.split('\t')[2]:
                    gene_name = line.split('\t')[8]
                    gene_name = gene_name.split(';')[0]
                    gene_name = gene_name.split('ID=')[1]
                    gene_name = gene_name.replace('"','')
                    
                    value = int(line.split('\t')[9])
                    
                    if gene_name not in values_dict:
                        values_dict[gene_name] = {}
                        gene_list.append(gene_name)
                        
                    if sample_uid not in values_dict[gene_name]:
                        values_dict[gene_name][sample_uid] = 0
                        
                    values_dict[gene_name][sample_uid] += value
        
        tab_file.close()
        
        return(values_dict, gene_list)
            
def build_coverage_table(cmd_dict):
    
    outfile_name = args.outfile_name
    outfile = open(outfile_name, 'w')
    
    values_dict = {}
    gene_list = []
    sample_uid_list = []
    
    set_name = cmd_dict['set_name']
    output_dir = cmd_dict['output_dir']
    
    genome_feature_file = cmd_dict['genome_features']
    if (genome_feature_file[-3:] == 'gff') or ('gff' in genome_feature_file[-5:]):
        runmode = 'gff'
    if (genome_feature_file[-3:] == 'gtf') or ('gtf' in genome_feature_file[-5:]):
        runmode = 'gtf'
        
    for strain_py in cmd_dict['samples']:
        if args.use_replicates:
            for replicate_py in cmd_dict['samples'][strain_py]:
                sample_uid = ('{}-{}').format(strain_py, replicate_py)
                
                tab_file_name = ('{output_dir}/{set_name}-{strain_py}-{replicate_py}.coverage.tab').format(
                output_dir = output_dir,
                set_name = set_name,
                strain_py = strain_py, 
                replicate_py = replicate_py)
                
                if sample_uid not in sample_uid_list:
                    sample_uid_list.append(sample_uid)
                        
                values_dict, gene_list = parse_coverage_file(tab_file_name, values_dict, gene_list, sample_uid, runmode)
            
        if not args.use_replicates:
            sample_uid = ('{}').format(strain_py)
            
            tab_file_name = ('{output_dir}/{set_name}-{strain_py}.coverage.tab').format(
            output_dir = output_dir,
            set_name = set_name,
            strain_py = strain_py)
            
            if sample_uid not in sample_uid_list:
                sample_uid_list.append(sample_uid)
                        
            values_dict, gene_list = parse_coverage_file(tab_file_name, values_dict, gene_list, sample_uid, runmode)
    
    gene_list.sort()
    sample_uid_list.sort()    
    
    header = '#gene'
    
    for sample_uid in sample_uid_list:
        header+=('\t'+sample_uid)
        
    header = header + '\n'
    
    outfile.write(header)
    
    for gene_name in gene_list:
        outline = gene_name 
        for sample_uid in sample_uid_list:
            if (gene_name in values_dict):
                if sample_uid in values_dict[gene_name]:
                    outline += '\t' + str(values_dict[gene_name][sample_uid])
                else:
                    outline += '\t' + str(np.nan)
                
        outline = outline + '\n'
                
        outfile.write(outline)

    outfile.close()
    
# ___  Body ___
    
cmd_dict, process = load_cmd_file()

if process and args.align:
    make_text(cmd_dict)
    
if process and args.rarefaction:
    run_rarefaction(cmd_dict)

if args.evaluate:
    summarize_stats(cmd_dict)
    
if args.coverage:
    coverage_run(cmd_dict)
    
if args.table_coverage_build:
    build_coverage_table(cmd_dict)
    
    
                    
    