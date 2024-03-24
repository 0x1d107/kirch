gnuplot << END
set yrange [*:*] reverse
set term png
set output "$2"
plot '$1' matrix with image;
END
