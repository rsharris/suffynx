# discordant mates coverage depth signal track. Dependency: perl, samtools, bedtools

(1) get the >=40 MAPQ and >=40% mapped (anchoring):

```bash  
    for which in distant broken bichromosomal ; do
        time bedtools bamtobed -cigar -i alignments/ZEB_150_150_MP.${which}.name_sorted.bam \
          | perl -lane \
            'print if $F[4] >=40' \
          | perl -lane \
            '@C = split ('\''([a-zA-Z])'\'', $F[6]);
             while (@C)
               {
               $l=shift(@C);
               $s=shift(@C);
               $clip=$clip+$l if $s=~m/[SH]$/;
               $mat=$mat+$l if $s=~m/[MDX=]$/
               }
             if ($clip+$mat)
               {print if $mat/($clip+$mat) >= 0.4}
             $clip=0;
             $mat=0' \
          > temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.bed
        done
```

(2) Extract their unmapped mates and apply the BBB and ### filter:

```bash  
    for which in distant broken bichromosomal ; do
        time (cat temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.bed;
              samtools view alignments/ZEB_150_150_MP.${which}.name_sorted.bam) \
          | perl -lane \
            'unless ($F[9])
               {
               @G=split/\//,$F[3];
               $H{$G[0]}=$G[1];
               $J{$G[0]}=$_
               }
             if ($F[9])
               {
               @E=split/\t/,$J{$F[0]} ;
               $mate=1 if $F[1] & 64;
               $mate=2 if $F[1] & 128;
               print "\@$F[0]\t$mate\n$F[9]\n\+\n$F[10]"
                 if $H{$F[0]}
                 and $H{$F[0]} ne $mate
                 and $F[10]!~m/^B{3,}/
                 and $F[10]!~m/B{3,}$/
                 and $F[10]!~m/^\#{3,}/
                 and $F[10]!~m/\#{3,}$/}' \
          > temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.fq
        done
```

(3) Q20 trim, BBB and ### trim and get rid of the short ones:

```bash  
    for which in distant broken bichromosomal ; do
        time seqtk trimfq -q 0.01 temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.fq \
          | seqtk seq -L 30 - \
          > temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.MSQF.fq
        done
```

(4a) Find the location of the trimmed anchoring mates. First generate a temporary sam using the mapping infor in the original bam file and the trimmed sequence in the fq file.

```bash  
    for which in distant broken bichromosomal ; do
        time (cat temp.ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.MSQF.fq;
              samtools view alignments/ZEB_150_150_MP.${which}.name_sorted.bam) \
          | perl -lane \
            'unless ($F[9])
               {
               $F[0]=~s/^\@//;
               $N{$F[0]}=$F[1] if $. % 4 eq 1
               }
             if ($F[9])
               {
               $p=1 if $F[1] & 64;
               $p=2 if $F[1] & 128;
               print if $N{$F[0]} and $N{$F[0]} ne $p
               }' \
          > temp.ZEB_150_150_MP.${which}.mate.sam.temp
        done
```

(4b) Remove duplicates, and then convert back to bed format.

```bash  
    for which in distant broken bichromosomal ; do
        time \
          (samtools view -h alignments/ZEB_150_150_MP.${which}.name_sorted.bam | grep ^@;
           cat temp.ZEB_150_150_MP.${which}.mate.sam.temp) \
          > temp.ZEB_150_150_MP.${which}.mate.bam
          samtools rmdup -S temp.ZEB_150_150_MP.${which}.mate.bam temp.ZEB_150_150_MP.${which}.mate.rmdup.bam
          | samtools view -Sbh temp.ZEB_150_150_MP.${which}.mate.rmdup.bam \
          | bedtools bamtobed -cigar -i - \
          > tracks/ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.MSQF.bed
        done
```

(5) sort and convert to bedgraph:

```bash  
    for which in distant broken bichromosomal ; do
        time bedtools sort -i tracks/ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.MSQF.bed \
          | bedtools genomecov -bg -g from_lei/hg19.chrlen.txt -i - \
          > tracks/ZEB_150_150_MP.${which}.name_sorted.MQ40.CP40.MSQF.bedgraph
        done
```

(6) combine them

```bash  
    trackName="ZEB_150_150_MP BDB MAPQ40 CLIP40 NoBBBPPP"
    echo "track type=bedGraph name=~${trackName}~ description=~${trackName}~" \
      | tr '~' '"' \
      > tracks/ZEB_150_150_MP.BDB.MQ40.CP40.MSQF.bedgraph
    time (cat tracks/ZEB_150_150_MP.broken.name_sorted.MQ40.CP40.MSQF.bed;
          cat tracks/ZEB_150_150_MP.distant.name_sorted.MQ40.CP40.MSQF.bed;
          cat tracks/ZEB_150_150_MP.bichromosomal.name_sorted.MQ40.CP40.MSQF.bed) \
      | bedtools sort -i - \
      | bedtools genomecov -bg -g \
          from_lei/hg19.chrlen.txt -i - \
      >> tracks/ZEB_150_150_MP.BDB.MQ40.CP40.MSQF.bedgraph
    ln -s        ZEB_150_150_MP.BDB.MQ40.CP40.MSQF.bedgraph \
          tracks/ZEB_150_150_MP.BDB.MQ40.CP40.MSQF.dat
```

(7) duplicate removal with samtools

```bash  
    for readtype in MP
      do
      for patient in S1 S2 S3 S4 S5 S6 S7 S8
        do
        for disctype in broken bichromosomal distant
          do
          samtools rmdup -S ${patient}_LGL_${readtype}.${disctype}.MMQ40.MCP40.bam ${patient}_LGL_${readtype}.${disctype}.MMQ40.MCP40.rmdup.bam
          done
        done
      done
```
    
(8) conversion to bedgraph

```bash  
    for readtype in MP
      do
      for patient in S1 S2 S3 S4 S5 S6 S7 S8
        do
        (for disctype in broken bichromosomal distant; do bedtools bamtobed -i ${patient}_LGL_${readtype}.${disctype}.MMQ40.MCP40.rmdup.bam; done) | bedtools sort -i - | bedtools genomecov -bg -g ~/leiyang/DB/hg19.chrlen.txt -i - > ${patient}_LGL_${readtype}.BDB.MMQ40.MCP40.rmdup.bedgraph
        done
      done
```
