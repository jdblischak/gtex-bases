#!/bin/bash

featureCounts \
  -a data/ENSG00000223745.saf \
  -o data/ENSG00000223745-chr1-93262186-93345888.txt \
  -F SAF \
  -f \
  --read2pos 5 \
  /project2/mstephens/dongyue/gtex/SRR597952_1_93297593_93307481.bam \
  /project2/mstephens/dongyue/gtex/SRR661409_20_57603732_57607422.bam
