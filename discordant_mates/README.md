# discordant mates coverage depth signal track. Dependency: perl, samtools, bedtools

(Lei working on this)

To fit in with the rest of the pipeline, this should
produce a track named tracks/ZEB_150_150_MP.rmdup.bedgraph

(1) get the >=40 MAPQ and >=40% mapped (anchoring):
    for which in distant broken bichromosomal ; do
        time bedtools bamtobed -cigar -i alignments/${run}.${which}.name_sorted.bam \
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
          > farf.${run}.${which}.name_sorted.MQ40.CP40.bed
        done

(2) Extract their unmapped mates and apply the BBB and ### filter:
    for which in distant broken bichromosomal ; do
        time (cat farf.${run}.${which}.name_sorted.MQ40.CP40.bed;
              samtools view alignments/${run}.${which}.name_sorted.bam) \
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
          > farf.${run}.${which}.name_sorted.MQ40.CP40.fq
        done

(3) Q20 trim, BBB and ### trim and get rid of the short ones:
    for which in distant broken bichromosomal ; do
        time seqtk trimfq -q 0.01 farf.${run}.${which}.name_sorted.MQ40.CP40.fq \
          | seqtk seq -L 30 - \
          > farf.${run}.${which}.name_sorted.MQ40.CP40.MSQF.fq
        done

(4a) Find the location of the trimmed anchoring mates. First generate a temporary sam using the mapping infor in the original bam file and the trimmed sequence in the fq file.
    for which in distant broken bichromosomal ; do
        time (cat farf.${run}.${which}.name_sorted.MQ40.CP40.MSQF.fq;
              samtools view alignments/${run}.${which}.name_sorted.bam) \
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
          > farf.${run}.${which}.mate.sam.temp
        done

(4b) Remove duplicates, and then convert back to bed format.
    today=`today {dd}{mmm}{yyyy}`
    for which in distant broken bichromosomal ; do
        time \
          (samtools view -h alignments/${run}.${which}.name_sorted.bam | grep ^@;
           cat farf.${run}.${which}.mate.sam.temp) \
          > farf.${run}.${which}.mate.bam
          samtools rmdup -S farf.${run}.${which}.mate.bam farf.${run}.${which}.mate.rmdup.bam
          | samtools view -Sbh farf.${run}.${which}.mate.rmdup.bam \
          | bedtools bamtobed -cigar -i - \
          > tracks/${run}.${which}.name_sorted.MQ40.CP40.MSQF.${today}.bed
        done

(5) sort and convert to bedgraph:
    today=`today {dd}{mmm}{yyyy}`
    for which in distant broken bichromosomal ; do
        time bedtools sort -i tracks/${run}.${which}.name_sorted.MQ40.CP40.MSQF.${today}.bed \
          | bedtools genomecov -bg -g from_lei/hg19.chrlen.txt -i - \
          > tracks/${run}.${which}.name_sorted.MQ40.CP40.MSQF.${today}.bedgraph
        done

(6) combine them
    today=`today {dd}{mmm}{yyyy}`
    trackName="${run} BDB MAPQ40 CLIP40 NoBBBPPP"
    echo "track type=bedGraph name=~${trackName}~ description=~${trackName}~" \
      | tr '~' '"' \
      > tracks/${run}.BDB.MQ40.CP40.MSQF.${today}.bedgraph
    time (cat tracks/${run}.broken.name_sorted.MQ40.CP40.MSQF.${today}.bed;
          cat tracks/${run}.distant.name_sorted.MQ40.CP40.MSQF.${today}.bed;
          cat tracks/${run}.bichromosomal.name_sorted.MQ40.CP40.MSQF.${today}.bed) \
      | bedtools sort -i - \
      | bedtools genomecov -bg -g \
          from_lei/hg19.chrlen.txt -i - \
      >> tracks/${run}.BDB.MQ40.CP40.MSQF.${today}.bedgraph
    ln -s        ${run}.BDB.MQ40.CP40.MSQF.${today}.bedgraph \
          tracks/${run}.BDB.MQ40.CP40.MSQF.${today}.dat

(7) duplicate removal with samtools
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
    
    
(8) conversion to bedgraph
    for readtype in MP
      do
      for patient in S1 S2 S3 S4 S5 S6 S7 S8
        do
        (for disctype in broken bichromosomal distant; do bedtools bamtobed -i ${patient}_LGL_${readtype}.${disctype}.MMQ40.MCP40.rmdup.bam; done) | bedtools sort -i - | bedtools genomecov -bg -g ~/leiyang/DB/hg19.chrlen.txt -i - > ${patient}_LGL_${readtype}.BDB.MMQ40.MCP40.rmdup.bedgraph
        done
      done

