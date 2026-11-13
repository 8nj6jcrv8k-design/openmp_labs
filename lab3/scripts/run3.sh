#!/bin/bash

SRC="lab3.f90"
N=300

echo "Компиляция версий -O0 и -O2..."
gfortran -O0 -fopenmp "$SRC" -o lab3_O0
gfortran -O2 -fopenmp "$SRC" -o lab3_O2

#Порядок вложенности циклов и опции оптимизации (-O0 и -O2)
echo "opt,method,threads,n,checksum,time" > res_p2_loops.csv
for OPT in O0 O2; do
    BIN="./lab3_${OPT}"
    for M in SERIAL serial_ijk serial_jki; do
        OUT=$(OMP_NUM_THREADS=1 $BIN $M $N)
        echo "${OPT},${OUT}" >> res_p2_loops.csv
    done
done

# Месторасположение директивы !$omp do (на 4 потоках, -O2)
echo "method,threads,n,checksum,time" > res_p5_placement.csv
for M in omp_outer omp_mid omp_inner; do
    OMP_NUM_THREADS=4 ./lab3_O2 $M $N >> res_p5_placement.csv
done

# Влияние OMP_SCHEDULE (static с разными блоками, dynamic, guided) на 4 потоках
echo "schedule,method,threads,n,checksum,time" > res_p6_schedule.csv
for SCH in "static" "static,1" "static,5" "static,20" "static,50" "dynamic" "guided"; do
    OUT=$(OMP_NUM_THREADS=4 OMP_SCHEDULE="$SCH" ./lab3_O2 omp_outer $N)
    echo "\"${SCH}\",${OUT}" >> res_p6_schedule.csv
done

#  Исследование collapse(1, 2, 3) при -O0 и -O2 и разных блоках static на 4 потоках
echo "opt,collapse_m,chunk,method,threads,n,checksum,time" > res_p7_collapse.csv
for OPT in O0 O2; do
    BIN="./lab3_${OPT}"
    for CHUNK in 1 10 100; do
        # m = 1 соответствует omp_outer
        OUT1=$(OMP_NUM_THREADS=4 OMP_SCHEDULE="static,${CHUNK}" $BIN omp_outer $N)
        echo "${OPT},1,${CHUNK},${OUT1}" >> res_p7_collapse.csv

        # m = 2
        OUT2=$(OMP_NUM_THREADS=4 OMP_SCHEDULE="static,${CHUNK}" $BIN collapse2 $N)
        echo "${OPT},2,${CHUNK},${OUT2}" >> res_p7_collapse.csv

        # m = 3
        OUT3=$(OMP_NUM_THREADS=4 OMP_SCHEDULE="static,${CHUNK}" $BIN collapse3 $N)
        echo "${OPT},3,${CHUNK},${OUT3}" >> res_p7_collapse.csv
    done
done

# Масштабируемость от 1 до 8 потоков (сравнение SERIAL, omp_outer и SECTIONS)
echo "method,threads,n,checksum,time" > res_p8_scaling.csv
OMP_NUM_THREADS=1 ./lab3_O2 SERIAL $N >> res_p8_scaling.csv
for T in 1 2 3 4 5 6 7 8; do
    OMP_NUM_THREADS=$T OMP_SCHEDULE="static" ./lab3_O2 omp_outer $N >> res_p8_scaling.csv
    OMP_NUM_THREADS=$T ./lab3_O2 SECTIONS $N >> res_p8_scaling.csv
done

echo "Все эксперименты завершены! Результаты сохранены в файлы res_p*.csv"
