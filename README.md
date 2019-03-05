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
```bash  
cd my_exe_path
ln -s some_path/suffynx/chrom_avg/close_intervals.py close_intervals
ln -s some_path/suffynx/chrom_avg/collect_tags.py collect_tags
ln -s some_path/suffynx/chrom_avg/create_script_call_insertions_discordant.py create_script_call_insertions_discordant
ln -s some_path/suffynx/chrom_avg/create_script_clipped_breakpoints.py create_script_clipped_breakpoints
ln -s some_path/suffynx/chrom_avg/create_script_clipped_breakpoints_high.py create_script_clipped_breakpoints_high
ln -s some_path/suffynx/chrom_avg/create_script_depth.py create_script_depth
ln -s some_path/suffynx/chrom_avg/create_script_discordant_mates_dense.py create_script_discordant_mates_dense
ln -s some_path/suffynx/chrom_avg/create_script_insert_depth.py create_script_insert_depth
ln -s some_path/suffynx/chrom_avg/create_script_insert_depth_dense.py create_script_insert_depth_dense
ln -s some_path/suffynx/chrom_avg/create_script_insert_depth_sparse.py create_script_insert_depth_sparse
ln -s some_path/suffynx/chrom_avg/create_script_insert_length.py create_script_insert_length
ln -s some_path/suffynx/chrom_avg/create_script_insert_length_sparse.py create_script_insert_length_sparse
ln -s some_path/suffynx/chrom_avg/create_script_insert_length_sparse_or_normal_inserts_sparse.py create_script_insert_length_sparse_or_normal_inserts_sparse
ln -s some_path/suffynx/chrom_avg/create_script_map.py create_script_map
ln -s some_path/suffynx/chrom_avg/create_script_short_or_discordant.py create_script_short_or_discordant
ln -s some_path/suffynx/chrom_avg/fill_genomic_interval_gaps.py fill_genomic_interval_gaps
ln -s some_path/suffynx/chrom_avg/filtered_sam_to_intervals.py filtered_sam_to_intervals
ln -s some_path/suffynx/chrom_avg/intervals_to_ucsc_catalog.py intervals_to_ucsc_catalog
ln -s some_path/suffynx/chrom_avg/keep_first.py keep_first
ln -s some_path/suffynx/chrom_avg/make_bigwig_info.py make_bigwig_info
ln -s some_path/suffynx/chrom_avg/make_bwa_jobs.py make_bwa_jobs
ln -s some_path/suffynx/chrom_avg/proximal_feature_intervals.py proximal_feature_intervals
ln -s some_path/suffynx/chrom_avg/sam_reader.py sam_reader
ln -s some_path/suffynx/chrom_avg/today.py today
```

# Pipeline Tutorial

We use a two-layer job paradigm. For a given sample (i.e. reads from an
individual) we first run a series of python programs to create job scripts.
Then we run the job scripts.

Generally each stage of the pipeline corresponds to a _track_, where a track
describes a value for each interval of the reference genome. There's a python
script for each stage, which will create a job script that will compute the
track for that stage. Job scripts are bash shell scripts.

The inputs to the process are a set of reads for the sample, a reference genome
and a pipeline control file. You can arrange these files any way you like, but
in the following example we assume that we have a single working directory with
the following subdirectories

* reads: This contains two fastq files for each sequencing run; usually this
will be one run for mate pair and another run for paired end.

* genomes: This contains the reference fasta file, the bwa index, a chromosome
lengths file, and any blacklist interval files.

* data: This contains the pipeline control file.

* jobs: The job scripts will be created here.

* alignments: Read-vs-reference alignments will be written here.

* tracks: Track files will be created here.

* temp: Temporary files will be created here.

In this example, the genome is named "reference", and the sample is named
"ZEB".
We start with paired end reads files reads/ZEB_PE.1.fastq and
reads/ZEB_PE.2.fastq,
and mate pair reads files reads/ZEB_MP.1.fastq and
reads/ZEB_MP.2.fastq.
data/control.dat has been copied from the
repository and modified if necessary.
There are two blacklist files, genomes/reference.Ns.dat and
genomes/repeat_masker.reference.dat.

## Create the scripts that will map the reads to the reference.

```bash  
create_script_map \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --ref="  {base}/genomes/reference.fa" \
      --reads="{base}/reads/{run}.{mate}.fastq" \
      --bam="  {base}/alignments/{run}" \
      --namesorted \
      --qualityfiltered \
  > jobs/ZEB_MP.map.sh
chmod +x jobs/ZEB_MP.map.sh
```

## Create the scripts that will compute the average mate pair insert length signal and indicator tracks.

This is called "Track 1" in supplementary methods step 2.

```bash  
create_script_insert_length \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
      --namesorted \
      --chroms={base}/genomes/reference.chrom_lengths \
      --track={base}/tracks/{run}.insert_length \
      --gzip \
  > jobs/ZEB_MP.insert_length.sh
chmod +x jobs/ZEB_MP.insert_length.sh
```

```bash  
create_script_insert_length_sparse \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --blacklist={base}/genomes/reference.Ns.dat \
      --blacklist={base}/genomes/repeat_masker.reference.dat \
      --input={base}/tracks/{run}.insert_length.gz \
      --track={base}/tracks/{run}.insert_length.sparse \
  > jobs/ZEB_MP.insert_length_sparse.sh
chmod +x jobs/ZEB_MP.insert_length_sparse.sh
```

## Create the scripts that will compute the mate pair short and normal insert coverage depth signal and indicator tracks.

These are called "Tracks 2 and 3" in supplementary methods step 3.

```bash  
create_script_insert_depth \
      --class=short \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
      --namesorted \
      --chroms={base}/genomes/reference.chrom_lengths \
      --track={base}/tracks/{run}.{kind}_inserts.depth \
  > jobs/ZEB_MP.insert_depth.short.sh
chmod +x jobs/ZEB_MP.insert_depth.short.sh
```

```bash  
create_script_insert_depth_dense \
      --class=short \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --input={base}/tracks/{run}.{kind}_inserts.depth \
      --track={base}/tracks/{run}.{kind}_inserts.depth.dense \
  > jobs/ZEB_MP.insert_depth_dense.sh
chmod +x jobs/ZEB_MP.insert_depth_dense.sh
```

```bash  
create_script_insert_depth \
      --class=normal \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
      --namesorted \
      --chroms={base}/genomes/reference.chrom_lengths \
      --track={base}/tracks/{run}.{kind}_inserts.depth \
  > jobs/ZEB_MP.insert_depth.normal.sh
chmod +x jobs/ZEB_MP.insert_depth.normal.sh
```

```bash  
create_script_insert_depth_sparse \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --blacklist={base}/genomes/reference.Ns.dat \
      --blacklist={base}/genomes/repeat_masker.reference.dat \
      --input={base}/tracks/{run}.{kind}_inserts.depth \
      --track={base}/tracks/{run}.{kind}_inserts.depth.sparse \
  > jobs/ZEB_MP.insert_depth_sparse.sh
chmod +x jobs/ZEB_MP.insert_depth_sparse.sh
```

## Create the scripts that will compute the mate pair discordant mates coverage depth signal and indicator tracks.

This is called "Track 4" in supplementary methods step 4.

See https://github.com/rsharris/suffynx/tree/master/discordant_mates for
creation of the signal track.

This creates the script that converts the signal track to an indicator track.

```bash  
create_script_discordant_mates_dense \
      --class=short \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --input={base}/tracks/{run}.BDB.MMQ40.MCP40.rmdup.bedgraph \
      --track={base}/tracks/{run}.discordant_mates.dense \
  > jobs/ZEB_MP.discordant_mates_dense.sh
chmod +x jobs/ZEB_MP.discordant_mates_dense.sh
```

## Create the scripts that will compute the paired end clipped breakpoints signal and indicator tracks.

This is called "Track 5" in supplementary methods step 5.

```bash  
create_script_clipped_breakpoints \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_PE \
      --bam={base}/alignments/{run}.ql_filtered.name_sorted.bam \
      --namesorted \
      --chroms={base}/genomes/reference.chrom_lengths \
      --track={base}/tracks/{run}.clipped_breakpoints \
  > jobs/ZEB_PE.clipped_breakpoints.sh
chmod +x jobs/ZEB_PE.clipped_breakpoints.sh
```

```bash  
create_script_clipped_breakpoints_high \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_PE \
      --chroms={base}/genomes/reference.chrom_lengths \
      --input={base}/tracks/{run}.clipped_breakpoints \
      --track={base}/tracks/{run}.clipped_breakpoints.high \
  > jobs/ZEB_PE.clipped_breakpoints_high.sh
chmod +x jobs/ZEB_PE.clipped_breakpoints_high.sh
```

## Create the scripts that will combine tracks 1 thru 5 and call insertions.

As per supplementary methods step 6.

```bash  
create_script_insert_length_sparse_or_normal_inserts_sparse \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --input1={base}/tracks/{run}.insert_length.sparse \
      --input2={base}/tracks/{run}.normal_inserts.depth.sparse \
      --track={base}/tracks/{run}.insert_length_sparse_or_normal_inserts_sparse \
  > jobs/ZEB_MP.track1_or_track2.sh
chmod +x jobs/ZEB_MP.track1_or_track2.sh
```

```bash  
create_script_short_or_discordant \
      --init=shebang:bash \
      --base="`pwd`" \
      ZEB_MP \
      --chroms={base}/genomes/reference.chrom_lengths \
      --input1={base}/tracks/{run}.short_inserts.depth.dense \
      --input2={base}/tracks/{run}.discordant_mates.dense \
      --track={base}/tracks/{run}.short_or_discordant \
  > jobs/ZEB_MP.track3_or_track4.sh
chmod +x jobs/ZEB_MP.track3_or_track4.sh
```

```bash  
create_script_call_insertions_discordant \
      --control=data/control.dat \
      --init=shebang:bash \
      --base="`pwd`" \
      peRun=ZEB_PE
      mpRun=ZEB_MP
      --chroms={base}/genomes/reference.chrom_lengths \
      --input="{base}/tracks/{mprun}.insert_length_sparse_or_normal_inserts_sparse" \
      --input="{base}/tracks/{mprun}.short_or_discordant" \
      --track="{base}/tracks/{run}.called_insertions" \
  > jobs/ZEB.called_insertions.sh
chmod +x jobs/ZEB.called_insertions.sh
```

## Compute.

Once all the jobs scripts have been created, they should be run, like this:

```bash  
./jobs/ZEB_MP.map.sh
./jobs/ZEB_MP.insert_length.sh
./jobs/ZEB_MP.insert_length_sparse.sh
./jobs/ZEB_MP.insert_depth.normal.sh
./jobs/ZEB_MP.insert_depth_sparse.sh
./jobs/ZEB_MP.insert_depth.short.sh
./jobs/ZEB_MP.insert_depth_dense.sh
./jobs/ZEB_MP.discordant_mates_dense.sh
./jobs/ZEB_PE.clipped_breakpoints.sh
./jobs/ZEB_PE.clipped_breakpoints_high.sh
./jobs/ZEB_MP.track1_or_track2.sh
./jobs/ZEB_MP.track3_or_track4.sh
./jobs/ZEB.called_insertions.sh
```


Feb/2019, Bob Harris (rsharris *at* bx *dot* psu *dot* edu)
