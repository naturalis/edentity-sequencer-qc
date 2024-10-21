# edentity-sequencer-qc

This repository contains shell and python scripts for quality control of sequencing data. The scripts are designed to 
be used to compare the quality of sequencing data from different sequencing platforms. We compare the following vendors:

- Element Bio AVITI - downloaded from https://nextnuc.gbiomed.kuleuven.be/index.php/s/mkLgECw2odysgaL (password
  protected, credentials available on request) - 20241011_AV242402_4854_1-RawData-4854.tar, unpacked into 
  /data/luka.lenaroto/eDentity/sequencer_tender/elements on MaaS 37
- Illumina - downloaded from naturalis@sftp.nrcnvwa.nl via sftp (credentials available on request) - 
  /data/luka.lenaroto/eDentity/sequencer_tender/241008_VH00147_37_AAG3LFWM5 on MaaS 37

## Data

### Elements

Total sum of *.fastq.gz file sizes: 43497579058 bytes (43.5 GB). Note that these are compressed by the
vendor, and we don't know the compression ratio. Hence, this is merely a very rough estimate.

Markers in multiplexed data set according to Elements:

- Empty
- ITS-Blank
- ITSmix-Blank
- ITSmix-NL2
- ITSmix-SB
- ITS-NL2
- ITS-SB
- LSU-Blank
- LSU-NL2
- LSU-SB
- SSU-Blank
- SSU-NL2
- SSU-SB

### Illumina

Total sum of *.fastq.gz file sizes: 115430885055 bytes (115.4 GB). Note that these are compressed by the
vendor, and we don't know the compression ratio. Hence, this is merely a very rough estimate.

Markers in multiplexed data set according to Illumina:

- Empty
- ITS-Blank1
- ITS-Blank2
- ITS-Blank3
- ITSmix-Blank1
- ITSmix-Blank2
- ITSmix-Blank3
- ITSmix-NL2
- ITSmix-SB
- ITS-NL2
- ITS-SB
- LSU-Blank1
- LSU-Blank2
- LSU-Blank3
- LSU-NL2
- LSU-SB
- SSU-Blank1
- SSU-Blank2
- SSU-Blank3
- SSU-NL2
- SSU-SB
- Undetermined

## Results

Intermediate results are stored in the [results](results) folder.