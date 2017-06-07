#!/bin/bash
awk -F '\,' '{if (NR!=1){print $1}}' result.csv | sort -m > l1.txt
awk -F '\,' '{if (NR!=1){print $1}}' FilteredGalaxyZoo40.csv | sort -m > l2.txt
awk '{ h[$0] = ! h[$0] } END { for (k in h) if (h[k]) print k }' l1.txt l2.txt > l3.txt
python remake.py l3.txt FilteredGalaxyZoo40.csv


