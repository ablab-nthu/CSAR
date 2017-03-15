#!/usr/bin/php
<?php
print "CSAR v1.0.0\n";
$totalTimeStart = getMicrotime();
$meg = "Usage: php csar.php [option] -t <target_contigs.fna file> -r <reference.fna file> [-nuc, -pro]";
$meg_h = "Using '-h' for help";
$outPath = ".";
foreach($argv as $k => $each){
	if($k == 0) continue;
	if(strstr($each, "-h")){
		print "$meg\n";
print <<<END
Option:
  -t <string>   Target genome (i.e., draft genome to be scaffolded)
  
  -r <string>   Reference genome
  
  -nuc          Use NUCmer to identify conserved genetic markers between target and reference genomes
  
  -pro          Use PROmer to identify conserved genetic markers between target and reference genomes
  
  -o <string>   Output folder to contain all the files returned by running CSAR
  
  -h            Show help message
END;
		print "\n";
		exit();
	}
	else if($each == "-t"){
		$draft1 = $argv[$k+1];
		if(@substr(trim(file_get_contents($draft1)), 0, 1) != ">"){
			die("$meg\n$meg_h\n");
		}
	}
	else if($each == "-r"){
		$draft2 = $argv[$k+1];
		if(@substr(trim(file_get_contents($draft2)), 0, 1) != ">"){
			die("$meg\n$meg_h\n");
		}
	}
	else if($each == "-nuc"){
		$mum = "nuc";
	}
	else if($each == "-pro"){
		$mum = "pro";
	}
	else if($each == "-o"){
		$outPath = $argv[$k+1];
		if($outPath == "") $outPath = ".";
	}
	else if($each == "-time"){
		$getRunTime = 1;
	}

}

if(!@mkdir($outPath) && !is_dir($outPath)){
	die("Can't create the output directory.\n");
}

$mummerRunTime = 0;
if(isset($mum)){
	if(strstr($draft2, "/")){
		$draft2CurrentPath = substr($draft2, strrpos($draft2, "/")+1);
	}
	else{
		$draft2CurrentPath = $draft2;
	}
	$mummerStartTime = getMicrotime();
	system($mum."mer $draft1 $draft2 -p $outPath/$draft2CurrentPath.".$mum);
	system("delta-filter -r -q $outPath/$draft2CurrentPath.$mum.delta > $outPath/$draft2CurrentPath.$mum.filter");
	system("show-coords -l -d $outPath/$draft2CurrentPath.$mum.filter > $outPath/$draft2CurrentPath.$mum.coords");
	$mummerRunTime = round(getMicrotime() - $mummerStartTime, 3);
	$coords = "$outPath/$draft2CurrentPath.$mum.coords";
}
else{
	die("$meg\n$meg_h\n");
}

if(!file_exists($coords) || !isset($coords)){
	die("$meg\n$meg_h\n");
}

//================================================
// Deal with .coords file
//================================================
$arr = file($coords, FILE_IGNORE_NEW_LINES);
if($arr[1] == "PROMER"){
	$mum = "pro";
}
else{
	$mum = "nuc";
}

for($i = 0; $i < 5; $i++){
	array_shift($arr);
}
$lines = array();

$element = preg_split("/ +/", trim($arr[0]));
$num = count($element);
if(($mum == "nuc" && $num != 17) || ($mum == "pro" && $num != 19)){
	die("$meg\n$meg_h\n");
}

foreach($arr as $k => $line){
	$element = preg_split("/ +/", trim($line));
	list($contig1_tag, $contig2_tag) = preg_split("/\t/", $element[$num-1]);
	$lines[$contig1_tag][$k][0] = $element[0]; //start 1
	$lines[$contig1_tag][$k][1] = $element[1]; //end 1
	if($lines[$contig1_tag][$k][0] > $lines[$contig1_tag][$k][1]){
		swap($lines[$contig1_tag][$k][0], $lines[$contig1_tag][$k][1]);	
	}
	$lines[$contig1_tag][$k][3] = $element[3]; //start 2
	$lines[$contig1_tag][$k][4] = $element[4]; //end 2
	if($lines[$contig1_tag][$k][3] > $lines[$contig1_tag][$k][4]){
		swap($lines[$contig1_tag][$k][3], $lines[$contig1_tag][$k][4]);	
	}
	$lines[$contig1_tag][$k][6] = $element[6]; //len 1
	$lines[$contig1_tag][$k][7] = $element[7]; //len 2
	$lines[$contig1_tag][$k][$num-3] = $element[$num-3]; //frm 1
	$lines[$contig1_tag][$k][$num-2] = $element[$num-2]; //frm 2
	$lines[$contig1_tag][$k][$num-1] = $contig2_tag;
}

// According to contig 1
// Sort by contig 1
foreach($lines as $contig_tag => $each_contig){
	foreach($each_contig as $i => $each_marker1){
		foreach($each_contig as $j => $each_marker2){
			if($j <= $i) continue;
			if($lines[$contig_tag][$i][0] > $lines[$contig_tag][$j][0]){
				swap($lines[$contig_tag][$i], $lines[$contig_tag][$j]);
			}
			else if($lines[$contig_tag][$i][0] == $lines[$contig_tag][$j][0]){
				if($lines[$contig_tag][$i][1] > $lines[$contig_tag][$j][1]){
					swap($lines[$contig_tag][$i], $lines[$contig_tag][$j]);
				}
			}
		}
	}
}

// Number the markers of contig 1
$n = 0;
$draft1_marker = "";
foreach($lines as $contig_tag => $each_contig){
	$draft1_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$n++;
		$draft1_marker .= $n."\n";	
		$lines[$contig_tag][$i]["marker"] = $n;
	}
	$draft1_marker .= "\n";	
}

// Copy markers form contig 1 to contig 2
$lines2 = array();
$n = 0;
foreach($lines as $contig_tag => $each_contig){
	foreach($each_contig as $i => $each_marker1){
		$contig2_tag = $lines[$contig_tag][$i][$num-1];
		$lines2[$contig2_tag][$n][0] = $lines[$contig_tag][$i][3]; // start
		$lines2[$contig2_tag][$n][1] = $lines[$contig_tag][$i][4]; // end
		if($lines[$contig_tag][$i][$num-3] * $lines[$contig_tag][$i][$num-2] > 0){
			$lines2[$contig2_tag][$n]["marker"] = $lines[$contig_tag][$i]["marker"]; // marker
		}
		else{
			$lines2[$contig2_tag][$n]["marker"] = -$lines[$contig_tag][$i]["marker"]; // marker
		}
		$n++;
	}
}

// Sort by contig 2
foreach($lines2 as $contig_tag => $each_contig){
	foreach($each_contig as $i => $each_marker1){
		foreach($each_contig as $j => $each_marker2){
			if($j <= $i) continue;
			if($lines2[$contig_tag][$i][0] > $lines2[$contig_tag][$j][0]){
				swap($lines2[$contig_tag][$i], $lines2[$contig_tag][$j]);
			}
			else if($lines2[$contig_tag][$i][0] == $lines2[$contig_tag][$j][0]){
				if($lines2[$contig_tag][$i][1] > $lines2[$contig_tag][$j][1]){
					swap($lines2[$contig_tag][$i], $lines2[$contig_tag][$j]);
				}
			}
		}
	}
}


// Record the markers of contig 2
$draft2_marker = "";
foreach($lines2 as $contig_tag => $each_contig){
	$draft2_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$draft2_marker .= $lines2[$contig_tag][$i]["marker"]."\n";
	}
	$draft2_marker .= "\n";	
}

//file_put_contents("$coords"."_1.contigs", $draft1_marker); 
//file_put_contents("$coords"."_2.contigs", $draft2_marker); 

//======================================
// Load markers of contigs
//======================================
$contigSet = array();
$contigSet[0] = explode(">", $draft1_marker);
$contigSet[1] = explode(">", $draft2_marker);
$originMarkers = $contigSet[0];

// Check markers
foreach($contigSet as $k => $set){
	$contigNum[$k] = 0;
	$total[$k] = 0;
	$allMarkers[$k] = array();
	$contigs[$k] = array();
	$contigTag[$k] = array();
	$contigBegin[$k] = array();
	$contigEnd[$k] = array();
	$telomere[$k] = array();
	$telomere_label[$k] = array();

	foreach($set as $l => $contig){
		if($contig == "") continue;
		$contigNum[$k]++;
		$markers = explode("\n", trim($contig));
		$markersNum = count($markers);
		$contigTag[$k][$markers[1]] = $markers[0];
		$contigTag[$k][$markers[$markersNum-1]] = $markers[0];
		array_push($contigBegin[$k], $markers[1]);
		array_push($contigEnd[$k], $markers[$markersNum-1] * -1); // head of a marker, e. g. "+3" for +3 (tail) and -3 (head) 
		array_push($telomere[$k], $markers[1]);
		array_push($telomere[$k], $markers[$markersNum-1] * -1);
		$telomere_label[$k][$markers[1]] = $l;
		$telomere_label[$k][$markers[$markersNum-1] * -1] = $l;

		array_shift($markers);
		$markers_tmp = $markers;
		array_walk($markers_tmp, "array_abs");
		$allMarkers[$k] = array_merge($allMarkers[$k], $markers_tmp);

		$total[$k] += count($markers);
		array_push($contigs[$k], $markers);
	}
	sort($allMarkers[$k]);
}

if($allMarkers[0] != $allMarkers[1]){
	die("The markers of two contig sets are not the same!\n");
}

foreach($allMarkers as $k => $markersSet){
	foreach($markersSet as $i => $marker){
		if(isset($markersSet[$i+1]) && $marker + 1 != $markersSet[$i+1]){
			die("There is a missing marker or repeat after the marker \"".$marker."\" of the input file ".($k+1)."!\n");
		}
	}
}

$contigSet = $contigs;

$CSAR_StartTime = getMicrotime();
//======================================
// 0. Add the complement strands
//======================================
foreach($contigSet as $k => $set){
	foreach($set as $l => $contig){
		$markersNum = count($contig);
		for($i = 0; $i < $markersNum; $i++){
			array_push($contigSet[$k][$l], $contig[$markersNum-1 - $i] * -1); 
		}
	}
}

//======================================
// 1. Get sigma-pi inverse
//======================================
foreach($contigSet[0] as $l => $contig){
	$reverse_contig[$l] = array_reverse($contig);
}

$sigma_piInverse = product($contigSet[1], $reverse_contig); // sigma and pi inverse
$sigma_piInverse_ori = $sigma_piInverse;

$index_in_sigma_piI = array();
foreach($sigma_piInverse_ori as $key => $cycle){
	if(count($cycle) < 2) continue;
	foreach($cycle as $marker){
		$index_in_sigma_piI[$marker] = $key;
	}
}

//======================================
// 2. Find good fusions to apply on pi and sigma
//======================================
// pi
foreach($telomere[0] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; // skip the 1-cycle 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum[0]*2; $j++){ // $contigNum[$k] * 2 telomeres
		if(!isset($telomere[0][$j])) continue;
		$y = $telomere[0][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; // skip the 1-cycle 
		// if there are x and y in c such that x and y belong to telomeres of pi and (x, y) can not devide pi
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[0][$x] != $telomere_label[0][$y]){
			$contigSet[0] = _join(array($x, $y), $contigSet[0]);
			$sigma_piInverse = _split($sigma_piInverse, array($x, $y));
			unset($telomere[0][$i]); // set the telomeres to null
			unset($telomere[0][$j]); // set the telomeres to null
			telomere_relabel($telomere_label[0], $x, $y);
		}
	}
}

// sigma
foreach($telomere[1] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; // skip the 1-cycle 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum[1]*2; $j++){ // $contigNum[$k] * 2 telomeres
		if(!isset($telomere[1][$j])) continue;
		$y = $telomere[1][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; // skip the 1-cycle 
		// if there are x and y in c such that x and y belong to telomeres of sigma and (x, y) can not devide sigma
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[1][$x] != $telomere_label[1][$y]){
			$contigSet[1] = _join(array($x, $y), $contigSet[1]);
			$sigma_piInverse = _split(array($x, $y), $sigma_piInverse);
			unset($telomere[1][$i]); // set the telomeres to null
			unset($telomere[1][$j]); // set the telomeres to null
			telomere_relabel($telomere_label[1], $x, $y);
		}
	}
}

$CSAR_RunTime = round(getMicrotime() - $CSAR_StartTime, 3);
//======================================
// 3. Output
//======================================
// Get the plus strand

$outPutTargetOnly = 1;

$plus_strand = array();
foreach($contigSet as $i => $eachSet){
	foreach($eachSet as $j => $eachContig){
		$markerNum = count($eachContig) / 2;
		foreach($eachContig as $k => $marker){
			foreach($telomere[$i] as $eachTelomere){
				if($marker == $eachTelomere){
					$plus_strand[$i][$j] = rotate_toBegin($eachContig, $k);
					array_splice($plus_strand[$i][$j], 0, $markerNum);
					break(2);
				}
			}
		}
	}
	if($outPutTargetOnly == 1){
		break;
	}
}

$plus_o = "0";
$minus_o = "1";
$mergeContig_tmp[0] = "";
$mergeContig_tmp[1] = "";

foreach($plus_strand as $i => $eachSet){
	$temp = "";
	foreach($eachSet as $j => $eachContig){
		$mergeContig_tmp[$i] .= ">Scaffold_".($j+1)."\n";
		foreach($eachContig as $marker){
			@$tagTmp = $contigTag[$i][$marker];
			@$tagTmpN = $contigTag[$i][$marker*-1];
			if(isset($tagTmp) && $tagTmp != $temp){
				$mergeContig_tmp[$i] .= $tagTmp." $plus_o\n";
				$temp = $tagTmp;
			}
			else if(isset($tagTmpN) && $tagTmpN != $temp){
				$mergeContig_tmp[$i] .= $tagTmpN." $minus_o\n";
				$temp = $tagTmpN;
			}
		}
		$mergeContig_tmp[$i] .= "\n";
	}
	$scaffoldNum = $j+1;
	
	if($i == 0){
		$outputFile = "scaffolds.$mum.csar";
	}
	else{
		$outputFile = "reference_scaffolds.$mum.csar";
	}
	file_put_contents("$outPath/$outputFile", $mergeContig_tmp[$i]);

	if($outPutTargetOnly == 1){
		break;
	}
}

// Generate the fasta seqence of scaffolds
$genSeqTime = getMicrotime();

$contigSeq = array();
parse_contig($draft1);
$list = file("$outPath/$outputFile");
$fasta = "";

foreach($list as $k => $each){
	if($each == "\n" || (strstr($each, ">") && ($fasta .= $each))){
		continue;
	}
	
	$partEach = substr($each, 0, strpos($each," ")+1);
	$partEach = trim($partEach);
	if(trim(substr($each, -2)) == $plus_o){
		$rewrite = $contigSeq[$partEach];
		$fasta .= $rewrite;
	}
	else{
		$rewrite = $contigSeq["$partEach.reverse"];
		$fasta .= $rewrite;
	}

	if($list[$k+1] == "\n" || strstr($list[$k+1], ">")){
		$fasta .= "\n";
	}
	else{
		$fasta .= str_repeat("N", 100)."\n";
	}

}
file_put_contents("$outPath/$outputFile.fna", $fasta);

$unmappedFile = 1;
if(isset($unmappedFile)){
	$unmappedContigs = "";
	$unmappedContigs_seq = "";
	$draftContig = file_get_contents($draft1);
	$draftContig = explode(">", $draftContig);
	$n = 0;
	foreach($draftContig as $eachInDraft){
		$toBeTok = $eachInDraft;
		$contigNameInD = strtok($toBeTok, " \n");
		foreach($originMarkers as $eachInContig){
			$contigNameInC = strtok($eachInContig, " \n");
			if($contigNameInD == $contigNameInC){
				$unmapped = 0;
				break;
			}
			else{
				$unmapped = 1;
			}
		}
		if($unmapped == 1){
			$scaffoldNum++;	
			$unmappedContigs .= ">Scaffold_".$scaffoldNum." (unmapped contig)\n";
			$unmappedContigs .= "$contigNameInD $plus_o\n\n";
			$unmappedContigs_seq .= ">Scaffold_".$scaffoldNum." (unmapped contig)\n";
			$unmappedContigs_seq .= $contigSeq[$contigNameInD]."\n";
		}
	}
	file_put_contents("$outPath/$outputFile", $unmappedContigs, FILE_APPEND);
	file_put_contents("$outPath/$outputFile.fna", $unmappedContigs_seq, FILE_APPEND);
}
$genSeqTime = round((getMicrotime() - $genSeqTime), 3);
$totalTime = round((getMicrotime() - $totalTimeStart), 3);
if(isset($getRunTime)){
	$runtimeContext = "Runtime info:\n";
	$runtimeContext .= "MUMmer: $mummerRunTime seconds\n";
	$runtimeContext .= "Main algorithm: $CSAR_RunTime seconds\n";
	$runtimeContext .= "Generate the fasta file: $genSeqTime seconds\n";
	$runtimeContext .= "Total Time: $totalTime seconds\n";
//	file_put_contents("csar_".$mum."_runtime.txt", $runtimeContext);
	print $runtimeContext;
}

//======================================
// Functions
//======================================
function array_abs(&$item, $key){
	$item = abs($item);
}

function swap(&$a, &$b){
	$temp = $a;
	$a = $b;
	$b = $temp;
}

function _split($a, $b){
	global $printWhile;

	if($printWhile == 1){
		print "in split\n";
		print_r($a);
		print_r($b);
	}

	if(!is_array($a[0])){
		foreach($b as $n => $each){
			$k1 = array_search($a[0], $each);
			if($k1 === false) continue;

			$new_arr = rotate_toBegin($each, $k1);
			$k2 = array_search($a[1], $new_arr);
			$part1 = array_slice($new_arr, 0, $k2);
			$part2 = array_slice($new_arr, $k2);
			array_splice($b, $n, 1, array($part1, $part2));
			break;
		}
	}
	else{
		$temp = $a;
		$a = $b;
		$b = $temp;
		foreach($b as $n => $each){
			$k1 = array_search($a[0], $each);
			if($k1 === false) continue;

			$new_arr = rotate_toEnd($each, $k1);
			$k2 = array_search($a[1], $new_arr);
			$part1 = array_slice($new_arr, 0, $k2+1);
			$part2 = array_slice($new_arr, $k2+1);
			array_splice($b, $n, 1, array($part1, $part2));
			break;
		}

	}

	if($printWhile == 1){
		print_r($b);
	}
	return $b;
}

function _join($a, $b){
	global $printWhile;

	if($printWhile == 1){
		print "in join\n";
		print_r($a);
		print_r($b);
	}

	if(!is_array($a[0])){
		foreach($b as $n => $each){
			$temp_k1 = array_search($a[0], $each);
			if($temp_k1 !== false){
				$n1 = $n;
				$k1 = $temp_k1;
			}
			$temp_k2 = array_search($a[1], $each);
			if($temp_k2 !== false){
				$n2 = $n;
				$k2 = $temp_k2;
			}

			if(isset($n1) && isset($n2)){
				$new_arr1 = rotate_toBegin($b[$n1], $k1);
				$new_arr2 = rotate_toBegin($b[$n2], $k2);
				$join_arr = array_merge($new_arr1, $new_arr2);
				break;
			}
		}
		$join_arr = array($join_arr);
		foreach($b as $n => $each){
			if($n == $n1 || $n == $n2) continue;
			array_push($join_arr, $each);
		}
	}
	else{
		$temp = $a;
		$a = $b;
		$b = $temp;
		foreach($b as $n => $each){
			$temp_k1 = array_search($a[0], $each);
			if($temp_k1 !== false || array($a[0]) == $each){
				$n1 = $n;
				$k1 = $temp_k1;
			}
			$temp_k2 = array_search($a[1], $each);
			if($temp_k2 !== false || array($a[1]) == $each){
				$n2 = $n;
				$k2 = $temp_k2;
			}

			if(isset($n1) && isset($n2)){
				$new_arr1 = rotate_toEnd($b[$n1], $k1);
				$new_arr2 = rotate_toEnd($b[$n2], $k2);
				$join_arr = array_merge($new_arr1, $new_arr2);
				break;
			}
		}
		$join_arr = array($join_arr);
		foreach($b as $n => $each){
			if($n == $n1 || $n == $n2) continue;
			array_push($join_arr, $each); 		
		}

	}
	if($printWhile == 1){
		print_r($join_arr);
	}
	return $join_arr;
}

function rotate_toBegin($b, $k){
	$n = count($b);
	$new_arr = array();
	for($i=0;$i<$n;$i++){
		if($k+$i >= $n){
			$new_arr[$i] = $b[$k+$i-$n];
		}
		else{
			$new_arr[$i] = $b[$k+$i];
		}
	}
	return $new_arr;
}

function rotate_toEnd($b, $k){
	$n = count($b);
	$new_arr = array();
	for($i=0;$i<$n;$i++){
		if($k+$i+1 >= $n){
			$new_arr[$i] = $b[$k+$i+1-$n];
		}
		else{
			$new_arr[$i] = $b[$k+$i+1];
		}
	}
	return $new_arr;
}

function product($a, $b){
	global $printWhile;

	if($printWhile == 1){	
		print "in product a:\n";
		print_r($a);
		print "in product b:\n";
		print_r($b);
	}

	$toSelf = array();
	$temp = array();
	foreach($b as $i => $m){
		$toSelf[$i] = $m;
		// Map to itself
		array_shift($toSelf[$i]);
		array_push($toSelf[$i],$m[0]);
		foreach($toSelf[$i] as $j => $n){
			foreach($a as $k => $o){
				$key = array_search($n, $o);
				if($key !== FALSE){
					if(isset($o[$key+1])){
						$temp[$i][$j] = $o[$key+1];
					}
					else{
						$temp[$i][$j] = $o[0]; // to the begin, cuz the circularity.
					}
				}
			}
		}
	}
	$merge_b = merge($b);
	$merge_temp = merge($temp);

	// Generate the product
	$prod = array();
	$temp_product = array();
	$i = 0;
	$k = 0;
	$flag = array();
	while($i<count($merge_b)){
		if(isset($flag[$i]) && $flag[$i] == 1){
			if($temp_product != array()){
				$prod[$k] = $temp_product;
				$temp_product = array();
				$k++;
			}
			$i++;
			continue;
		}
		array_push($temp_product, $merge_temp[$i]);
		$flag[$i] = 1;
		$i = array_search($merge_temp[$i], $merge_b);
	}
	if($printWhile == 1){
		print "product:\n";
		print_r($prod);
	}
	return $prod;
}

function merge($array){
	$merge = array();
	foreach($array as $i){
		$merge = array_merge($merge, $i);
	}
	return $merge;
}

function parse_contig($multiFasta){
	global $contigSeq;
	$lines = file($multiFasta, FILE_IGNORE_NEW_LINES);
	$forward = "";
	$reverse = "";
	$forward_arr = array();
	$reverse_arr = array();
	$pass = 0;

	foreach($lines as $k => $line){
		if($line == "") continue;
		
		if(substr($line, 0, 1) != ">"){
			array_push($forward_arr, $line."\n");
			array_push($reverse_arr, strrev($line)."\n");
		}

		if(substr($line, 0, 1) == ">" || !isset($lines[$k+1])){
			if($pass == 1){
				foreach($forward_arr as $each){
					$forward .= $each;
				}
				$n = count($reverse_arr);
				for($i = $n-1; $i>=0; $i--){
					$reverse .= $reverse_arr[$i];
				}

				$contigSeq[$fileName] = $forward;		
				$new_reverse = $reverse;
				$seqLen = strlen($reverse);
				for($i = 0; $i < $seqLen; $i++){
					if($reverse[$i] == "A"){
						$new_reverse[$i] = "T";
					}
					else if($reverse[$i] == "T"){
						$new_reverse[$i] = "A";
					}
					else if($reverse[$i] == "C"){
						$new_reverse[$i] = "G";
					}
					else if($reverse[$i] == "G"){
						$new_reverse[$i] = "C";
					}
				}
				$contigSeq["$fileName.reverse"] = $new_reverse;		
				$forward_arr = array();
				$reverse_arr = array();
				$forward = "";
				$reverse = "";
			}
			
			if(strstr($line, " ")){
				$fileName = substr($line, 1, strpos($line, " ")-1);
			}
			else{
				$fileName = substr($line, 1);
			}
			$pass = 1;
		}
	}
}

function telomere_relabel(&$telomere_label, $x, $y){
	$telomere_label_to_change = $telomere_label[$y];
	foreach($telomere_label as $key => $telomere_end){
		if($telomere_end == $telomere_label_to_change){
			$telomere_label[$key] = $telomere_label[$x];
		}
	}

}

function getMicrotime(){
	list($usec, $sec) = explode(" ", microtime());
	return ((float)$usec + (float)$sec);
}
?>
