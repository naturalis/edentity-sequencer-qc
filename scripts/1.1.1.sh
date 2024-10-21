# https://docs.google.com/spreadsheets/d/1qg99bWHMhyynXRqUq34QofuUcpc2gbsQ/edit?gid=714344309#gid=714344309
# Calculate the fraction of bases with quality score ≥ Q30 in a FASTQ file
FILES=$(find . -name "*.fastq.gz")
for FILE in $FILES; do
  PHRED=$(zcat $FILE | awk 'NR%4==0 {total += length; for(i=1; i<=length; i++) if(substr($0,i,1)>="?") q30++} END {print total, ",", q30/total}')
  echo $FILE,$PHRED
done
