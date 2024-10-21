# edentity-sequencer-qc

This repository contains shell and python scripts for quality control of sequencing data. The scripts are designed to 
be used to compare the quality of sequencing data from different sequencing platforms. We compare the following vendors:

- Element Bio AVITI - downloaded from https://nextnuc.gbiomed.kuleuven.be/index.php/s/mkLgECw2odysgaL (password
  protected, credentials available on request) - 20241011_AV242402_4854_1-RawData-4854.tar, unpacked into 
  /data/luka.lenaroto/eDentity/sequencer_tender/elements on MaaS 37
- Illumina - downloaded from naturalis@sftp.nrcnvwa.nl via sftp (credentials available on request) - 
  /data/luka.lenaroto/eDentity/sequencer_tender/241008_VH00147_37_AAG3LFWM5 on MaaS 37

Markers in multiplexed data set:

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
