<?php

$uR=0.0+$_GET['uR'];
$DR=0.0+$_GET['DR'];
$LifTR=0.0+$_GET['LifTR'];
$xrange=0.0+$_GET['xrange'];
$tmax=0.0+$_GET['tmax'];


header('Content-type: image/png');
echo `python3 mRNATimeEvolution.py $uR $DR $xrange $tmax $LifTR`;
?>