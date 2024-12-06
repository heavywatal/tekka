set -eux
TMP_OUT=tekka_cli.sh

./tekka -o $TMP_OUT
rm -r $TMP_OUT

./tekka --seed 42 -y40 -K2000 -r2 -l2 --sa 2,2 --sj 2,2 -o $TMP_OUT
rm -r $TMP_OUT
