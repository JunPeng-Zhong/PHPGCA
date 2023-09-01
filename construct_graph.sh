#!/bin/bash

RED="\033[31m"
GREEN="\033[32m"
SUFFIX="\033[0m"

echo -e  $RED"================== BUILD PROTEIN =================="$SUFFIX
python build_protein.py --addition-test True


echo -e  $RED"================== VIRUS-VIRUS =================="$SUFFIX
python build_contigs.py
python edge_virus_virus.py

echo -e  $RED"================== VIRUS-PROKARYOTE =================="$SUFFIX
python build_pairlink.py --refresh True --refresh-allHOST True
python build_edge_virus_prokaryote.py

echo -e  $RED"================== HETERO GRAPH =================="$SUFFIX
python create_feature_withoutfeat.py
python multimodal_graph2.py

echo -e  $GREEN"================== SUCCESS =================="$SUFFIX