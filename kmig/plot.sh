INDEX=${3:-1}
gnuplot << END
set yrange [*:*] reverse
set term png
set output "$2"
unset key
plot '$1' matrix with image;
END
