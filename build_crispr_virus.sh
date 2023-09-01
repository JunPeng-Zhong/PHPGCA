#!/bin/bash

echo "Build allCRISPRs & allVIRUS"

echo "Building crispr..."
python build_allCRISPR.py

echo "Building virus..."
python build_allVIRUS.py