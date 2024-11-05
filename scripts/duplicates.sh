samtools collate -@ 24 -O -u *phix-mapped-only.bam | \
samtools fixmate -@ 24 -m -u - - | \
samtools sort -@ 24 -u - | \
samtools markdup -@ 24 - markdup.bam \ 
samtools flagstat -@ 24 markdup.bam