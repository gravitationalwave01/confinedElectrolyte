set term pngcairo enhanced 
set output "profile.png"
set yrange[0:]
set title "Ion profile"
set xlabel "Distance / l_b"
set ylabel "Concentration {/Symbol r}*l_b^3"

plot "output.dat" u 1:2 w l title "cation",\
"output.dat" u 1:3 w l title "anion"


set output "pot.png"
set title "Electric potential"
set autoscale y
set xlabel "Distance / l_b"
set ylabel "Electric potential e{/Symbol b}{/Symbol y}"
plot "output.dat" u 1:4 w l notitle 
