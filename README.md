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
```

# Tutorial

(to be written soon)


Feb/2019, Bob Harris (rsharris *at* bx *dot* psu *dot* edu)
