#!/bin/bash

cat result_*.csv > result.csv
rm result_*.csv
rm temp_*_upstream.tsv
rm temp_*_downstream.tsv
