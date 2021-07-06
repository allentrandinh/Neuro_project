<?php

$uR=0.0+$_GET['uR'];
$DR=0.0+$_GET['DR'];
$LifTR=0.0+$_GET['LifTR'];
$DP=0.0+$_GET['DP'];
$uP=0.0+$_GET['uP'];
$LifTP=0.0+$_GET['LifTP'];
$xrange=0.0+$_GET['xrange'];
$tmax=0.0+$_GET['tmax'];
$SomDentr=0.0+$_GET['SomDentr'];


header('Content-type: image/png');
echo `python3 ProteinTimeEvolution_web.py $uR $DR $LifTR $SomDentr $DP $uP $LifTP $xrange $tmax`;
?>
