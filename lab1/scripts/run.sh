#!/usr/bin/env bash
# Использование: ./scripts/run.sh <МЕТОД> <выходной_файл.csv>
# Пример:        ./scripts/run.sh REDUCTION results/reduction.csv

set -e   # остановить скрипт при первой ошибке

# --- 1. Проверяем аргументы ---
if [ $# -ne 2 ]; then
    echo "Использование: $0 <МЕТОД> <выходной_файл.csv>" >&2
    exit 1
fi

METHOD=$1
OUT=$2

# --- 2. Настройки ---
# Путь к нашей скомпилированной программе прямо на Рабочем столе
EXE="./integral_fortran"
N=100000000            # число отрезков разбиения
THREADS="1 2 4 8"      # потоки, которые проверяем на Маке

# Последовательному методу потоки не нужны: достаточно одного прогона
if [ "$METHOD" = "SERIAL" ]; then
    THREADS="1"
fi

# --- 3. Проверяем, что программа собрана ---
if [ ! -x "$EXE" ]; then
    echo "Программа $EXE не найдена. Сначала соберите проект." >&2
    exit 1
fi

# --- 4. Основная часть: заголовок CSV и цикл по числу потоков ---
echo "method,threads,n,integral,time" > "$OUT"

for t in $THREADS; do
    echo "Запуск $METHOD на $t потоках..."
    OMP_NUM_THREADS=$t "$EXE" "$METHOD" "$N" >> "$OUT"
done

echo "Готово: результаты записаны в $OUT"
