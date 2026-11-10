#!/usr/bin/env bash
set -e   

if [ $# -ne 2 ]; then
    echo "Использование: $0 <МЕТОД> <выходной_файл.csv>" >&2
    exit 1
fi

METHOD=$1
OUT=$2
EXE="./vector_program"
N=100000000            # Рекомендуемый размер массива (100 млн элементов)
THREADS="1 2 4 8"      

if [ "$METHOD" = "SERIAL" ]; then
    THREADS="1"
fi

if [ ! -x "$EXE" ]; then
    echo "Программа $EXE не найдена. Сначала скомпилируйте код." >&2
    exit 1
fi

echo "method,threads,n,vector_sum,time" > "$OUT"

for t in $THREADS; do
    echo "Запуск $METHOD на $t потоках..."
    OMP_NUM_THREADS=$t "$EXE" "$METHOD" "$N" >> "$OUT"
done

echo "Готово: результаты в $OUT"
