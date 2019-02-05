# suffynx
Support scripts for "Retrovirus insertion site analysis of LGL leukemia patient
genomes"

# Dependencies

* genodsp, which can be found at https://github.com/rsharris/genodsp

# Installation

To install suffynx from source code:  
1. Download the latest version of suffynx using Github  
```bash  
    cd some_path  
    git clone https://github.com/rsharris/suffynx  
```  
2. Compile:  
```bash  
    cd some_path/suffynx/chrom_avg  
    make  

    cd some_path/suffynx/encodachrom  
    make  
```

3. Set up paths:  

Copy or symlink the executables (some_path/suffynx/chrom_avg/chrom_avg and
some_path/suffynx/encodachrom/encodachrom) somewhere into your shell's $PATH.
For example
```bash  
    cd my_exe_path
    ln -s some_path/suffynx/chrom_avg/chrom_avg .
    ln -s some_path/suffynx/encodachrom/encodachrom .
```

For every python script, symlink it somewhere into your shell's $PATH, without
the ".py" extensions. For example
```bash  
    cd my_exe_path
    ln -s some_path/suffynx/chrom_avg/close_intervals.py close_intervals
    ln -s some_path/suffynx/chrom_avg/collect_tags.py collect_tags
    ln -s some_path/suffynx/chrom_avg/create_script_call_insertions_discordant.py create_script_call_insertions_discordant
    ln -s some_path/suffynx/chrom_avg/create_script_depth.py create_script_depth
    ln -s some_path/suffynx/chrom_avg/create_script_discordant_mates_dense.py create_script_discordant_mates_dense
    ln -s some_path/suffynx/chrom_avg/create_script_insert_depth.py create_script_insert_depth
    ln -s some_path/suffynx/chrom_avg/create_script_insert_depth_dense.py create_script_insert_depth_dense
    ln -s some_path/suffynx/chrom_avg/create_script_insert_depth_sparse.py create_script_insert_depth_sparse
    ln -s some_path/suffynx/chrom_avg/create_script_insert_length.py create_script_insert_length
    ln -s some_path/suffynx/chrom_avg/create_script_insert_length_sparse.py create_script_insert_length_sparse
    ln -s some_path/suffynx/chrom_avg/create_script_insert_length_sparse_or_normal_inserts_sparse.py create_script_insert_length_sparse_or_normal_inserts_sparse
    ln -s some_path/suffynx/chrom_avg/create_script_short_or_discordant.py create_script_short_or_discordant
    ln -s some_path/suffynx/chrom_avg/fill_genomic_interval_gaps.py fill_genomic_interval_gaps
    ln -s some_path/suffynx/chrom_avg/filtered_sam_to_intervals.py filtered_sam_to_intervals
    ln -s some_path/suffynx/chrom_avg/intervals_to_ucsc_catalog.py intervals_to_ucsc_catalog
    ln -s some_path/suffynx/chrom_avg/keep_first.py keep_first
    ln -s some_path/suffynx/chrom_avg/make_bigwig_info.py make_bigwig_info
    ln -s some_path/suffynx/chrom_avg/make_bwa_jobs.py make_bwa_jobs
    ln -s some_path/suffynx/chrom_avg/proximal_feature_intervals.py proximal_feature_intervals
    ln -s some_path/suffynx/chrom_avg/sam_reader.py sam_reader
```

# Pipeline Tutorial

(to be completed soon)

……… inputs
……… reads, alignments
……… directory layout alignments, data, tracks, jobs
……… copy data/control.dat

Create the script that will compute the average mate pair insert length signal
track.

```bash  
    create_script_insert_length \
          --control=data/control.dat \
          --init=shebang:bash \
          --base="`pwd`" \
          SAMPLE_ID \
          --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
          --namesorted \
          --chroms={base}/data/hg19.chrom_lengths \
          --track={base}/tracks/{run}.insert_length \
          --gzip \
      > jobs/SAMPLE_ID.insert_length.sh
```

Feb/2019, Bob Harris (rsharris *at* bx *dot* psu *dot* edu)
