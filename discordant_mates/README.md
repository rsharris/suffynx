## Discordant mates coverage depth signal track.

Dependency: perl, samtools, bedtools, seqtk.

(1) Filter alignments:

```bash  
samtools view -h alignments/ZEB_MP.name_sorted.bam \
  | filtered_sam_to_intervals --namesorted \
      --requiremates \
      --prohibit:"(CIGAR == *)" \
      --prohibit:"(RNEXT == =)" \
      --require:" (FLAGS & 0xF00 == 0)" \
      --justsamrecords \
  | samtools view -bS - \
  > alignments/ZEB_MP.bichromosomal.name_sorted.bam

samtools view -h alignments/ZEB_MP.name_sorted.bam \
  | filtered_sam_to_intervals --namesorted \
      --requiremates \
      --require:"broken mates with unmapped" \
      --justsamrecords \
  | samtools view -bS - \
  > alignments/ZEB_MP.broken.name_sorted.bam

samtools view -h alignments/ZEB_MP.name_sorted.bam \
  | filtered_sam_to_intervals --namesorted \
      --requiremates=30K.. \
      --prohibit:"(CIGAR == *)" \
      --require:" (RNEXT == =)" \
      --require:" (FLAGS & 0xF00 == 0)" \
      --justsamrecords \
  | samtools view -bS - \
  > alignments/ZEB_MP.distant.name_sorted.bam
```

(2) Get the >=40 MAPQ and >=40% mapped (anchoring):

```bash  
for disctype in broken bichromosomal distant ; do
    bedtools bamtobed -cigar -i alignments/ZEB_MP.${disctype}.name_sorted.bam \
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
      > temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.bed
    done
```

(3) Extract their unmapped mates and apply the BBB and ### filter:

```bash  
for disctype in broken bichromosomal distant ; do
    (cat temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.bed;
     samtools view alignments/ZEB_MP.${disctype}.name_sorted.bam) \
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
      > temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.fq
    done
```

(4) Q20 trim, BBB and ### trim and get rid of the short ones:

```bash  
for disctype in broken bichromosomal distant ; do
    seqtk trimfq -q 0.01 temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.fq \
      | seqtk seq -L 30 - \
      > temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.MSQF.fq
    done
```

(5a) Find the location of the trimmed anchoring mates. First generate a
temporary sam using the mapping info in the original bam file and the trimmed
sequence in the fq file.

```bash  
for disctype in broken bichromosomal distant ; do
    (cat temp.ZEB_MP.${disctype}.name_sorted.MQ40.CP40.MSQF.fq;
     samtools view alignments/ZEB_MP.${disctype}.name_sorted.bam) \
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
      > temp.ZEB_MP.${disctype}.mate.sam
    done
```

(5b) Remove duplicates, and then convert back to bed format.

```bash  
for disctype in broken bichromosomal distant ; do
    (samtools view -h alignments/ZEB_MP.${disctype}.name_sorted.bam | grep ^@;
     cat temp.ZEB_MP.${disctype}.mate.sam) \
      > temp.ZEB_MP.${disctype}.mate.bam
    samtools rmdup \
        -S temp.ZEB_MP.${disctype}.mate.bam \
        temp.ZEB_MP.${disctype}.mate.rmdup.bam
      | samtools view -Sbh temp.ZEB_MP.${disctype}.mate.rmdup.bam \
      | bedtools bamtobed -cigar -i - \
      > tracks/ZEB_MP.${disctype}.name_sorted.MQ40.CP40.MSQF.bed
    done
```

(6) sort and convert to bedgraph:

```bash  
for disctype in broken bichromosomal distant ; do
    bedtools sort -i tracks/ZEB_MP.${disctype}.name_sorted.MQ40.CP40.MSQF.bed \
      | bedtools genomecov -bg -g hg19.chrlen.txt -i - \
      > tracks/ZEB_MP.${disctype}.name_sorted.MQ40.CP40.MSQF.bedgraph
    done
```

(7) combine them

```bash  
trackName="ZEB_MP BDB MAPQ40 CLIP40 NoBBBPPP"
echo "track type=bedGraph name=~${trackName}~ description=~${trackName}~" \
  | tr '~' '"' \
  > tracks/ZEB_MP.BDB.MQ40.CP40.MSQF.bedgraph

(cat tracks/ZEB_MP.broken.name_sorted.MQ40.CP40.MSQF.bed;
 cat tracks/ZEB_MP.distant.name_sorted.MQ40.CP40.MSQF.bed;
 cat tracks/ZEB_MP.bichromosomal.name_sorted.MQ40.CP40.MSQF.bed) \
  | bedtools sort -i - \
  | bedtools genomecov -bg -g hg19.chrlen.txt -i - \
  >> tracks/ZEB_MP.BDB.MQ40.CP40.MSQF.bedgraph
```

(8) duplicate removal with samtools

```bash  
for disctype in broken bichromosomal distant ; do
    samtools rmdup -S ZEB_MP.${disctype}.MMQ40.MCP40.bam \
                      ZEB_MP.${disctype}.MMQ40.MCP40.rmdup.bam
    done
```

(9) conversion to bedgraph

```bash  
for disctype in broken bichromosomal distant; do \
    bedtools bamtobed -i ZEB_MP.${disctype}.MMQ40.MCP40.rmdup.bam
    done \
  | bedtools sort -i - \
  | bedtools genomecov -bg -g hg19.chrlen.txt -i - \
  > ZEB_MP.BDB.MMQ40.MCP40.rmdup.bedgraph
```
