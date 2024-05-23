#!/bin/bash

for i in '0001' '0002' '0003' '0004' '0005' '0006' '0007' '0008' '0009' '0010'
do
  
  echo "Running data set $i"
	unzip -p track1a_20220404.zip patient_year/acic_patient_year_$i.csv > patient_year_$i.csv
	unzip -p track1a_20220404.zip patient/acic_patient_$i.csv > patient_$i.csv
	unzip -p track1a_20220404.zip practice_year/acic_practice_year_$i.csv > practice_year_$i.csv
	unzip -p track1a_20220404.zip practice/acic_practice_$i.csv > practice_$i.csv

	Rscript acic_mixed_trimming $i >> results_$i.txt

	rm patient_year_$i.csv
	rm patient_$i.csv
	rm practice_year_$i.csv
	rm practice_$i.csv
done



