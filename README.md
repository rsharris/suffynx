# suffynx
Support scripts for "Retrovirus insertion site analysis of LGL leukemia patient
genomes"

# Dependencies

* genodsp, which can be found at https://github.com/rsharris/genodsp
* bwa (we used version 0.7.12-r1039)

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
(NEED TO UPDATE THIS LIST)
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

We use a two-layer job paradigm. For a given sample (i.e. reads from an
individual) we first run a series of python programs to create job scripts.
Then we run the job scripts.

Generally each stage of the pipeline corresponds to a _track_, where a track
describes a value for each interval of the reference genome. There's a python
script for each stage, which will create a job script that will compute the
track for that stage.  Job scripts are bash shell scripts.

The inputs to the process are a set of reads for the sample, a reference genome
and a pipeline control file. You can arrange these files any way you like, but
in the following example we assume that we have a single working directory with
the following subdirectories

* reads: This contains two fastq files for each sequencing run; usually this
will be one run for paired end and another run for mate pair.

* genomes: This contains the reference fasta file, the bwa index, a chromosome
lengths file, and any blacklist interval files.

* data: This contains the pipeline control file.

* jobs: The job scripts will be created here.

* alignments: Read-vs-reference alignments will be written here.

* tracks: Track files will be created here.

* temp: Temporary files will be created here.

In this example, the genome is named "reference", and the sample is named
"ZEB".

……… copy control.dat from the repo



```bash  
    :> data/pe_run_names
    echo "SIM_101_101_PE" >> data/pe_run_names

    :> data/mp_run_names
    echo "SIM_150_150_MP" >> data/mp_run_names
```

Create the scripts that will map the reads to the reference.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_map \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --ref="  {base}/genomes/reference.fa" \
                --reads="{base}/reads/{run}.infected_reads.{mate}.fastq" \
                --bam="  {base}/alignments/{run}" \
                --namesorted \
                --qualityfiltered \
            > jobs/${run}.map.sh
          chmod +x jobs/${run}.map.sh
          done
```

Create the scripts that will compute the average mate pair insert length signal
tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_length \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
                --namesorted \
                --chroms={base}/genomes/reference.chrom_lengths \
                --track={base}/tracks/{run}.insert_length \
                --gzip \
            > jobs/${run}.insert_length.sh
          chmod +x jobs/${run}.insert_length.sh
          done
```

Create the scripts that will compute the average mate pair insert length
indicator tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_length_sparse \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --chroms={base}/genomes/reference.chrom_lengths \
                --blacklist={base}/genomes/reference.Ns.dat \
                --blacklist={base}/genomes/repeat_masker.reference.dat \
                --input={base}/tracks/{run}.insert_length.gz \
                --track={base}/tracks/{run}.insert_length.sparse \
            > jobs/${run}.insert_length_sparse.sh
          chmod +x jobs/${run}.insert_length_sparse.sh
          done
```

Create the scripts that will compute the normal insert coverage depth signal
tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_depth \
                --class=normal \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
                --namesorted \
                --chroms={base}/genomes/reference.chrom_lengths \
                --track={base}/tracks/{run}.{kind}_inserts.depth \
            > jobs/${run}.insert_depth.normal.sh
          chmod +x jobs/${run}.insert_depth.normal.sh
          done
```

Create the scripts that will compute the normal insert coverage depth
indicator tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_depth_sparse \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --chroms={base}/genomes/reference.chrom_lengths \
                --blacklist={base}/genomes/reference.Ns.dat \
                --blacklist={base}/genomes/repeat_masker.reference.dat \
                --input={base}/tracks/{run}.{kind}_inserts.depth \
                --track={base}/tracks/{run}.{kind}_inserts.depth.sparse \
            > jobs/${run}.insert_depth_sparse.sh
          chmod +x jobs/${run}.insert_depth_sparse.sh
          done
```

Create the scripts that will compute the short insert coverage depth signal
tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_depth \
                --class=short \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
                --namesorted \
                --chroms={base}/genomes/reference.chrom_lengths \
                --track={base}/tracks/{run}.{kind}_inserts.depth \
            > jobs/${run}.insert_depth.short.sh
          chmod +x jobs/${run}.insert_depth.short.sh
          done
```

Create the scripts that will compute the short insert coverage depth
indicator tracks.

```bash  
    cat data/pe_run_names data/mp_run_names \
      | while read run ; do
          create_script_insert_depth_dense \
                --class=short \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --chroms={base}/genomes/reference.chrom_lengths \
                --input={base}/tracks/{run}.{kind}_inserts.depth \
                --track={base}/tracks/{run}.{kind}_inserts.depth.dense \
            > jobs/${run}.insert_depth_dense.sh
          chmod +x jobs/${run}.insert_depth_dense.sh
          done
```

Create the scripts that will compute the mate pair discordant mates coverage
depth signal tracks.

```bash  
???
uses create_script_discordant_mates_dense on just the mp_run_names,
but there's some pre-processing that has to happen first
???
```

Create the scripts that will compute the mate pair discordant mates coverage
depth indicator tracks.

```bash  
(Bob to add this)
```

Create the scripts that will compute the paired end clipped breakpoints signal
tracks.

```bash  
    cat data/pe_run_names data/pe_run_names \
      | while read run ; do
          create_script_clipped_breakpoints \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
                --namesorted \
                --chroms={base}/genomes/reference.chrom_lengths \
                --track={base}/tracks/{run}.clipped_breakpoints \
            > jobs/${run}.clipped_breakpoints.sh
          chmod +x jobs/${run}.clipped_breakpoints.sh
          done
```

Create the scripts that will compute the paired end clipped breakpoints
indicator tracks.

```bash  
    cat data/pe_run_names data/pe_run_names \
      | while read run ; do
          create_script_clipped_breakpoints_high \
                --control=data/control.dat \
                --init=shebang:bash \
                --base="`pwd`" \
                ${run} \
                --chroms={base}/genomes/reference.chrom_lengths \
                --input={base}/tracks/{run}.clipped_breakpoints \
                --track={base}/tracks/{run}.clipped_breakpoints.high \
            > jobs/${run}.clipped_breakpoints_high.sh
          chmod +x jobs/${run}.clipped_breakpoints_high.sh
          done
```

(Bob to add the combination track stuff)

Feb/2019, Bob Harris (rsharris *at* bx *dot* psu *dot* edu)
