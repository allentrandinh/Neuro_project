<?php

$DP=0.0+$_GET['DP'];

$uP=0.0+$_GET['uP'];

$LifTP=0.0+$_GET['LifTP'];

$Dur=0.0+$_GET['Dur'];

$xrange=0.0+$_GET['xrange'];

$tmax=0.0+$_GET['tmax'];

$x0=0.0+$_GET['x0'];
header('Content-type: image/png');

echo `python3 Function3_plot_only.py $DP $uP $tmax $xrange $Dur $x0 $LifTP`;

?>

