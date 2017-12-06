#!/usr/bin/env php
<?php
/* CSAR v1.1.1 */
ini_set('memory_limit', -1);
$totalTimeStart = getMicrotime();
$opt = getopt("t:r:o:h", array('nuc', 'pro', 'dotplot', 'time', 'coords:', 'c1:', 'c2:', 'debug', 'more'));
$meg = 'Usage: php csar.php [option] -t <target_contigs.fna file> -r <reference.fna file> [--nuc, --pro]';
$meg_h = 'Using option \'-h\' for help';

if(isset($opt['h'])){
	print <<<END
$meg
Option:
	-t <string>   Target genome (i.e., draft genome to be scaffolded)

	-r <string>   Reference genome

	--nuc         Use NUCmer to identify conserved genetic markers between target and reference genomes

	--pro         Use PROmer to identify conserved genetic markers between target and reference genomes

	-o <string>   Output folder to contain all the files returned by running CSAR (default: ./csar_out)

	-h            Show help message

END;
	exit();
}

$CWD = getcwd();
$target = isset($opt['t']) ? realpath($opt['t']) : '';
$ref = isset($opt['r']) ? realpath($opt['r']) : '';
if($target != '') checkFASTA($target, 't');
if($ref != '') checkFASTA($ref, 'r');

$mum = isset($opt['nuc']) ? 'nuc' : '';
if($mum == '') $mum = isset($opt['pro']) ? 'pro' : '';
if(isset($opt['nuc']) && isset($opt['pro']))
	pErr('Either NUCmer or PROmer but not both should be used at once.');

$dotplot = isset($opt['dotplot']) ? 1 : 0;
$getRunTime = isset($opt['time']) ? 1 : 0;

$coords = isset($opt['coords']) ? realpath($opt['coords']) : '';
if($coords != '')
	$useCoords = substr($coords, -6) == 'coords' ? 1 : pErr($meg_h);

$inputContig1 = isset($opt['c1']) ? realpath($opt['c1']) : '';
if($inputContig1 != '')
	$useContig = substr($inputContig1, -7) == 'contigs' ? 1 : pErr($meg_h);

$inputContig2 = isset($opt['c2']) ? realpath($opt['c2']) : '';
if($inputContig2 != '')
	$useContig = substr($inputContig2, -7) == 'contigs' ? 1 : pErr($meg_h);

$outPath = isset($opt['o']) ? $opt['o'] : $CWD.'/csar_out';

if(!@mkdir($outPath, 0777, true) && !is_dir($outPath)){
	pErr('Can not create the output directory.');
}

$debug = isset($opt['debug']) ? 1 : 0;
$moreFile = isset($opt['more']) ? 1 : 0;

if(isset($useContig)){
	$mummerRunTime = 0;
	$mumOption = '';
	$mumOptSuf = 'test';
	echo "useContig\n";
	goto useContig;
}	
elseif(isset($useCoords)){
	$mummerRunTime = 0;
	$mumOption = "";
}
elseif($mum != ''){
	$refCurrentPath = $outPath.'/'.basename($ref);
	
	$mumOption = '';
	$mumOptSuf = $mum;
	$deltaPara = '-r -q'; // -r: one-to-many; -q: many-to-one
	$mummerRunTime = 0;
	$mummerStartTime = getMicrotime();
	system($mum."mer $mumOption $ref $target -p $refCurrentPath.".$mumOptSuf);
	system("delta-filter $deltaPara $refCurrentPath.$mumOptSuf.delta > $refCurrentPath.$mumOptSuf.filter"); 
	system("show-coords -l -d $refCurrentPath.$mumOptSuf.filter > $refCurrentPath.$mumOptSuf.coords");
	$coords = "$refCurrentPath.$mumOptSuf.coords";
	$mummerRunTime = round(getMicrotime() - $mummerStartTime, 3);
}
else{
	pErr($meg_h);
}

if(!file_exists($coords) || $coords == ''){
	pErr($meg_h);
}
//================================================
// Deal with .coords file
//================================================
$arr = file($coords, FILE_IGNORE_NEW_LINES);
$mum = ($arr[1] == 'PROMER') ? 'pro' : 'nuc';
$mumOptSuf = $mum;

for($i = 0; $i < 5; $i++){
	array_shift($arr);
}
$lines = array();

$element = preg_split('/ +/', trim($arr[0]));
$num = count($element);
if(($mum == 'nuc' && $num != 17) || ($mum == 'pro' && $num != 19)){
	pErr($meg_h);
}

$weight = 0;
$getMarkerTime = array();
$getMarkerTime[0] = getMicrotime();
foreach($arr as $k => $line){
	$element = preg_split('/ +/', trim($line));

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
$getMarkerTime[0] = round(getMicrotime() - $getMarkerTime[0], 3);
$getMarkerTime[1] = getMicrotime();

// According to contig 1
// Sort by contig 1
foreach($lines as $contig_tag => $each_contig){
	usort($lines[$contig_tag], 'cmp');
}

$getMarkerTime[1] = round(getMicrotime() - $getMarkerTime[1], 3);
$getMarkerTime[2] = getMicrotime();
// Number the markers of contig 1
$n = 0;
$ref_marker = '';
foreach($lines as $contig_tag => $each_contig){
	$ref_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$n++;
		$ref_marker .= $n."\n";	
		$lines[$contig_tag][$i]['marker'] = $n;
	}
	$ref_marker .= "\n";	
}
$getMarkerTime[2] = round(getMicrotime() - $getMarkerTime[2], 3);

$getMarkerTime[3] = getMicrotime();
// Copy markers form contig 1 to contig 2
$lines2 = array();
$n = 0;
foreach($lines as $contig_tag => $each_contig){
	foreach($each_contig as $i => $each_marker1){
		$contig2_tag = $lines[$contig_tag][$i][$num-1];
		$lines2[$contig2_tag][$n][0] = $lines[$contig_tag][$i][3]; // start
		$lines2[$contig2_tag][$n][1] = $lines[$contig_tag][$i][4]; // end
		if($lines[$contig_tag][$i][$num-3] * $lines[$contig_tag][$i][$num-2] > 0){
			$lines2[$contig2_tag][$n]['marker'] = $lines[$contig_tag][$i]['marker']; // marker
		}
		else{
			$lines2[$contig2_tag][$n]['marker'] = -$lines[$contig_tag][$i]['marker']; // reversed marker
		}
		$n++;
	}
}
$getMarkerTime[3] = round(getMicrotime() - $getMarkerTime[3], 3);

$getMarkerTime[4] = getMicrotime();
// Sort by contig 2
foreach($lines2 as $contig_tag => $each_contig){
	usort($lines2[$contig_tag], 'cmp');
}
$getMarkerTime[4] = round(getMicrotime() - $getMarkerTime[4], 3);

$getMarkerTime[5] = getMicrotime();
// Record the markers of contig 2
$target_marker = '';
foreach($lines2 as $contig_tag => $each_contig){
	$target_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$target_marker .= $lines2[$contig_tag][$i]['marker']."\n";
	}
	$target_marker .= "\n";	
}
$getMarkerTime[5] = round(getMicrotime() - $getMarkerTime[5], 3);

$coords = basename($coords);
if($moreFile == 1){
	file_put_contents("$outPath/$coords".'_ref.contigs', $ref_marker); 
	file_put_contents("$outPath/$coords".'_tar.contigs', $target_marker); 
}

//p("getMarkerTime:");
//print_r($getMarkerTime);
useContig:
//======================================
// Load markers of contigs
//======================================
$contigSet = array();

if(isset($useContig)){
	$ref_marker = file_get_contents($inputContig1);
	$target_marker = file_get_contents($inputContig2);
}
$contigSet[0] = explode('>', $ref_marker);
$contigSet[1] = explode('>', $target_marker);
$ref_originMarkers = $contigSet[0];
$tar_originMarkers = $contigSet[1];

// Check markers
foreach($contigSet as $k => $set){
	$contigNum[$k] = 0;
	$contigs[$k] = array();
	$contigTag[$k] = array();
	$telomere[$k] = array();
	$telomere_label[$k] = array();

	$contig_idx = 0;
	foreach($set as $l => $contig){
		if($contig == '') continue;
		$contigNum[$k]++;
		$markers = explode("\n", trim($contig));
		$markersNum = count($markers);
		$contigTag[$k][$markers[1]] = $markers[0];
		$contigTag[$k][$markers[$markersNum-1]] = $markers[0];
		array_push($telomere[$k], $markers[1]);
		array_push($telomere[$k], $markers[$markersNum-1] * -1);
		$telomere_label[$k][$markers[1]] = $contig_idx;
		$telomere_label[$k][$markers[$markersNum-1] * -1] = $contig_idx;
		$contig_idx++;
		array_shift($markers);
		array_push($contigs[$k], $markers);
	}
}

$contigSet = $contigs;
$stepTime = array();
$CSAR_StartTime = getMicrotime();
//======================================
// 0. Add the complement strands
//======================================
$stepTime[0] = getMicrotime();
foreach($contigSet as $k => $set){
	foreach($set as $l => $contig){
		$markersNum = count($contig);
		for($i = 0; $i < $markersNum; $i++){
			array_push($contigSet[$k][$l], $contig[$markersNum-1 - $i] * -1); 
		}
	}
}
$stepTime[0] = round(getMicrotime() - $stepTime[0], 3);
//p("stepTime 0: $stepTime[0]");
//======================================
// 1. Get sigma-pi inverse
//======================================
$stepTime[1] = getMicrotime();
foreach($contigSet[0] as $l => $contig){
	$reverse_contig[$l] = array_reverse($contig);
}

$sigma_piInverse = product($contigSet[1], $reverse_contig); // sigma and pi inverse

$index_in_sigma_piI = array();
foreach($sigma_piInverse as $key => $cycle){
	if(count($cycle) < 2) continue;
	foreach($cycle as $marker){
		$index_in_sigma_piI[$marker] = $key;
	}
}
$stepTime[1] = round(getMicrotime() - $stepTime[1], 3);
//p("stepTime 1: $stepTime[1]");

//======================================
// 2. Find good fusions to apply on pi and sigma
//======================================
// pi
$stepTime[2] = getMicrotime();
if(0)
foreach($telomere[0] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; // skip the 1-cycle 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum[0]*2; $j++){ // $contigNum[$k] * 2 telomeres
		if(!isset($telomere[0][$j])) continue;
		$y = $telomere[0][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; // skip the 1-cycle 
		// if there are x and y in c such that x and y belong to telomeres of pi and (x, y) can not devide pi
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[0][$x] != $telomere_label[0][$y]){
			$contigSet[0] = _join(array($x, $y), $contigSet[0], $telomere_label[0][$x], $telomere_label[0][$y]);
			unset($telomere[0][$i]); // set the telomeres to null
			unset($telomere[0][$j]); // set the telomeres to null
			telomere_relabel($telomere_label[0], $x, $y);
		}
	}
}

// sigma
$contigNum_temp = $contigNum[1] << 1; // $contigNum[$k] * 2 telomeres
foreach($telomere[1] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; // skip the 1-cycle 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum_temp; $j++){ 
		if(!isset($telomere[1][$j])) continue;
		$y = $telomere[1][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; // skip the 1-cycle 
		// if there are x and y in c such that x and y belong to telomeres of sigma and (x, y) can not devide sigma
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[1][$x] != $telomere_label[1][$y]){
			if($debug == 1){
				p("telomere_label 1:");
				print_r($telomere_label[1]);
			}
			$contigSet[1] = _join(array($x, $y), $contigSet[1], $telomere_label[1][$x], $telomere_label[1][$y]);
			unset($telomere[1][$i]); // set the telomeres to null
			unset($telomere[1][$j]); // set the telomeres to null
			telomere_relabel($telomere_label[1], $x, $y);
		}
	}
}

$contigSet[1] = array_values($contigSet[1]);

$stepTime[2] = round(getMicrotime() - $stepTime[2], 3);
//p("stepTime 2: $stepTime[2]");

$CSAR_RunTime = round(getMicrotime() - $CSAR_StartTime, 3);
//======================================
// 3. Output
//======================================
// Get the plus strands
$stepTime[3] = getMicrotime();
$plus_strand = array();
foreach($contigSet[1] as $j => $eachContig){
	$markerNum = count($eachContig) / 2;
	foreach($eachContig as $k => $marker){
		foreach($telomere[1] as $eachTelomere){
			if($marker == $eachTelomere){
				$plus_strand[1][$j] = rotate_toBegin($eachContig, $k);
				array_splice($plus_strand[1][$j], $markerNum);
				break(2);
			}
		}
	}
}

$plus_o = '0';
$minus_o = '1';
$scaffoldNum = 0;

$temp = '';
foreach($plus_strand[1] as $j => $eachContig){
	$k = $j+1;
	$major_ori = 0;
	$scaffold[$k] = array();
	foreach($eachContig as $marker){
		@$tagTmp = $contigTag[1][$marker];
		@$tagTmpN = $contigTag[1][$marker*-1];
		if(isset($tagTmp) && $tagTmp != $temp){
			array_push($scaffold[$k], "$tagTmp $plus_o");
			$temp = $tagTmp;
			$major_ori++;
		}
		else if(isset($tagTmpN) && $tagTmpN != $temp){
			array_push($scaffold[$k], "$tagTmpN $minus_o");
			$temp = $tagTmpN;
			$major_ori--;
		}
	}
	//flip scaffolds to be plus in majority
	if($major_ori < 0){
		$scaffold[$k] = flip_scaffold($scaffold[$k]);
	}
}
$scaffoldNum = $j+1;

usort($scaffold, 'cmp_arrSize');

$all_scaffolds = '';
foreach($scaffold as $i => $scaffold_each){
	$all_scaffolds .= '>Scaffold_'.($i+1)."\n";
	foreach($scaffold_each as $line){
		$all_scaffolds .= "$line\n";
	}
	$all_scaffolds .= "\n";
}

$outputFile_target = "$outPath/scaffolds.$mumOptSuf.csar";
file_put_contents($outputFile_target, $all_scaffolds);

$stepTime[3] = round(getMicrotime() - $stepTime[3], 3);
//p("stepTime 3: $stepTime[3]");

// Generate the fasta seqence of scaffolds
$genSeqTime = getMicrotime();

$contigSeq = array();
if($target != ''){
	parse_contig($target);

	$genFnaFile = 1;
	if($genFnaFile == 1){
		genFnaFile($target, $outputFile_target);
	}
	
	$genUnmapped = 1;
	if($genUnmapped == 1){
		genUnmappedFile($target, $tar_originMarkers, $scaffoldNum, $outputFile_target, '');
	}
	$genSeqTime = round((getMicrotime() - $genSeqTime), 3);

	if($dotplot == 1 && $ref != ''){
		$dotplotStartTime = getMicrotime();
		$assembFile = '';
		$draft_assembly = $outputFile_target;
		$assembled_fasta = "$draft_assembly.fna";

		system($mum."mer $mumOption $ref $assembled_fasta -p $outPath/plotInput");
		system("delta-filter $deltaPara $outPath/plotInput.delta > $outPath/plotInput.filter");
		system("show-coords -l -d $outPath/plotInput.filter > $outPath/plotInput.coords");
		$direct = judge_dir("$outPath/plotInput.coords");
		if($direct == '-1'){
			$list = file($draft_assembly);
			$list = array_reverse($list);
			$list_r = array();
			$fasta = '';
			foreach($list as $each){
				if(trim($each) == '') continue;
				if(strstr($each, '>')){
					if(!strstr($each, "\n")){
						$each .= "\n";
					}
					array_push($list_r, $each);
					continue;
				}
				if(trim(substr($each, -2)) == $minus_o){
					array_push($list_r, str_replace(" $minus_o\n", " $plus_o\n", $each));
				}
				else{
					array_push($list_r, str_replace(" $plus_o\n", " $minus_o\n", $each));
				}
			}

			$replace_posi = 0;
			$flag = 0;
			foreach($list_r as $k => $each){
				if($each[0] == '>'){
					if($replace_posi == 0){ // the first tag
						$first_tag = $each;
						if(!isset($list_r[$k+1])){ // only one scaffold
							array_unshift($list_r, $first_tag);
							array_pop($list_r);
							$flag = 1;
							break;
						}
						$replace_posi = $k;
						continue;
					}
					$list_r[$replace_posi] = "\n".$list_r[$k]; // replace
					$replace_posi = $k;
				}
			}
			if($flag == 0){
				array_unshift($list_r, $first_tag);
				array_pop($list_r);
			}
			foreach($list_r as $each){
				$assembFile .= $each;
			}
			$assembFile .= "\n";

			file_put_contents("$draft_assembly.reverse", $assembFile);

			genFnaFile($target, "$draft_assembly.reverse");

			system($mum."mer $mumOption $ref $draft_assembly.reverse.fna -p $outPath/plotInput");
		}

		system("mummerplot -f --large --png $outPath/plotInput.delta -p $draft_assembly".".dotplot_HR");
		$gpFile = file_get_contents($draft_assembly.".dotplot_HR.gp");
		$gp = str_replace(' tiny ', ' large ', $gpFile);
		$gp = str_replace('1400,1400', '1600,1400', $gp);
		$gp = str_replace("\"QRY\"\n", "\"Assembly\"\n", $gp);
		file_put_contents($draft_assembly.".dotplot_HR.gp", $gp);
		system("gnuplot $draft_assembly".".dotplot_HR.gp");
		$dotplotRunTime = round(getMicrotime() - $dotplotStartTime, 2);
	}
}

$totalTime = round((getMicrotime() - $totalTimeStart), 3);
if($getRunTime == 1){
	$runtimeContext = "Runtime info:\n";
	$runtimeContext .= "MUMmer: $mummerRunTime seconds\n";
	$runtimeContext .= "Main algorithm: $CSAR_RunTime seconds\n";
	$runtimeContext .= "Generate the fasta file: $genSeqTime seconds\n";
	$runtimeContext .= "Total Time: $totalTime seconds\n";
	p($runtimeContext);
	file_put_contents("$outPath/csar_$mum.runTime", $runtimeContext);
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

function _join($a, $b, $x_label, $y_label){
	global $debug;

	if($debug == 1){
		p("in join");
		print_r($a);
		print_r($b);
		p("x_label: $x_label");
		p("y_label: $y_label");
	}


	$k1 = array_search($a[0], $b[$x_label]);
	$k2 = array_search($a[1], $b[$y_label]);

//	print "k1: $k1\n";
//	print "k2: $k2\n";

	$new_arr1 = rotate_toBegin($b[$x_label], $k1);
	$new_arr2 = rotate_toBegin($b[$y_label], $k2);
	$join_arr = array_merge($new_arr1, $new_arr2);
	unset($b[$x_label]);
	unset($b[$y_label]);
	$b[$x_label] = $join_arr;

	if($debug == 1){
		p("after join");
		print_r($b);
	}
	return $b;
}

function rotate_toBegin($b, $k){
	return array_merge(array_slice($b, $k), array_slice($b, 0, $k));
}

function rotate_toEnd($b, $k){
	return array_merge(array_slice($b, $k+1), array_slice($b, 0, $k+1));
}

function product($a, $b){
	global $debug;

	if($debug == 1){
		p("in product a:");
		print_r($a);
		p("in product b:");
		print_r($b);
	}
	
//	$runTime = getMicrotime();

	$a_mapping = array();
	$b_mapping = array();
	$temp = array();
	foreach($a as $i => $markers){
		$temp = $markers;
		// map to itself
		array_shift($temp);
		array_push($temp, $markers[0]);
		$a_mapping += array_combine($markers, $temp);
	}
//	print "a_mapping\n";
//	print_r($a_mapping);

	$numMarkers = 0;
	$temp = array();
	foreach($b as $i => $markers){
		$temp = $markers;
		// map to itself
		array_shift($temp);
		array_push($temp,$markers[0]);
		foreach($temp as $j => $n){
			$b_mapping[$markers[$j]] = $a_mapping[$n];
			$numMarkers++;
		}
	}
//	print "ProductRunTime: ".(getMicrotime() - $runTime)."\n";

	// generate the product
	$prod = array();
	$temp = array();
	$i = 0;
//	p("numMarkers: $numMarkers");
	while($i < $numMarkers){
	//	print "i: $i\n";
		foreach($b_mapping as $k => $marker){
			$m = isset($nextMarker) ? $nextMarker : $k;
			array_push($temp, $m);
			$i++;
			$nextMarker = $b_mapping[$m];
			unset($b_mapping[$m]);
			if(!isset($b_mapping[$nextMarker]) && $temp != array()){
				array_push($prod, $temp);
				$temp = array();
				unset($nextMarker);
				break;
			}
		}
	}

	if($debug == 1){
		p("product:");
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
	$lines = array();
	$handle = fopen($multiFasta, 'r');
	while(($line = fgets($handle)) !== false){
		array_push($lines, rtrim($line));
	}
	fclose($handle);

	$forward = '';
	$reverse = '';
	$forward_arr = array();
	$reverse_arr = array();
	$pass = 0;

	foreach($lines as $k => $line){
		if($line == '') continue;
		
		if($line[0] != '>'){
			array_push($forward_arr, $line."\n");
			array_push($reverse_arr, strrev($line)."\n");
		}
		
		if($line[0] == '>' || !isset($lines[$k+1])){
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
					if($reverse[$i] == 'A'){
						$new_reverse[$i] = 'T';
					}
					else if($reverse[$i] == 'T'){
						$new_reverse[$i] = 'A';
					}
					else if($reverse[$i] == 'C'){
						$new_reverse[$i] = 'G';
					}
					else if($reverse[$i] == 'G'){
						$new_reverse[$i] = 'C';
					}
					else if($reverse[$i] == 'a'){
						$new_reverse[$i] = 't';
					}
					else if($reverse[$i] == 't'){
						$new_reverse[$i] = 'a';
					}
					else if($reverse[$i] == 'c'){
						$new_reverse[$i] = 'g';
					}
					else if($reverse[$i] == 'g'){
						$new_reverse[$i] = 'c';
					}
				}
				$contigSeq["$fileName.reverse"] = $new_reverse;
				$forward_arr = array();
				$reverse_arr = array();
				$forward = '';
				$reverse = '';
			}
			
			if(strstr($line, ' ')){
				$fileName = substr($line, 1, strpos($line, ' ')-1);
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
	$keys = array_keys($telomere_label, $telomere_label_to_change);
	foreach($keys as $key){
		$telomere_label[$key] = $telomere_label[$x];
	}

}

function genFnaFile($draftFileName, $outputFile_draft){
	global $contigSeq;
	global $plus_o;
	
	//$contig_source = "contig_source";
	//$parseContigTime = getMicrotime();
	//parse_contig($draftFileName);
	//$parseContigTime = round(getMicrotime() - $parseContigTime, 3);
	
	$list = file($outputFile_draft);
	$append = '';
	$edge = 0;
	$outSuffix = '.fna';
	file_put_contents($outputFile_draft.$outSuffix, '');
	foreach($list as $k => $each){
		if($each == "\n" || strstr($each, '>')){
			file_put_contents($outputFile_draft.$outSuffix, $each, FILE_APPEND);
			$edge = 1;
			continue;
		}

		if($edge == 0){
		//	$append = isset($list[$k+1]) && $list[$k+1] != "\n" ? str_repeat('N', 100)."\n" : "\n";
			$append = str_repeat('N', 100)."\n";
			file_put_contents($outputFile_draft.$outSuffix, $append, FILE_APPEND);
		}

		$partEach = trim(substr($each, 0, strpos($each,' ')+1));
		if(trim(substr($each, -2)) == $plus_o){
			file_put_contents($outputFile_draft.$outSuffix, $contigSeq[$partEach], FILE_APPEND);
		}
		else{
			file_put_contents($outputFile_draft.$outSuffix, $contigSeq["$partEach.reverse"], FILE_APPEND);
		}
		
		$edge = 0;
	}
}

function genUnmappedFile($draftFileName, $originMarkers, $scaffoldNum, $outputFile_draft, $refPrefix){
	global $plus_o;
	global $contigSeq;

	$unmappedContigs = '';
	$unmappedContigs_seq = '';

	$draftContig = array();
	$handle = fopen($draftFileName, 'r');
	while(($line = fgets($handle)) !== false){
		if($line[0] == '>'){
			array_push($draftContig, $line);
		}
	}
	fclose($handle);
	$n = 0;
	foreach($draftContig as $eachInDraft){
		$toBeTok = $eachInDraft;
		$contigNameInD = strtok(substr($toBeTok, 1), " \n");
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
		//	$unmappedContigs .= ">".$refPrefix."Scaffold_".$scaffoldNum."_(unmapped_contig)\n";
			$unmappedContigs .= ">".$refPrefix."Scaffold_".$scaffoldNum."\n";
			$unmappedContigs .= "$contigNameInD $plus_o\n\n";
		//	$unmappedContigs_seq .= ">".$refPrefix."Scaffold_".$scaffoldNum."_(unmapped_contig)\n";
			$unmappedContigs_seq .= ">".$refPrefix."Scaffold_".$scaffoldNum."\n";
			$unmappedContigs_seq .= $contigSeq[$contigNameInD]."\n";
		}
	}
	file_put_contents($outputFile_draft, $unmappedContigs, FILE_APPEND);
	file_put_contents("$outputFile_draft.fna", $unmappedContigs_seq, FILE_APPEND);
}

function genUnmappedFile_joinAll($draftFileName, $originMarkers, $scaffoldNum, $outputFile_draft, $refPrefix){
	global $plus_o;
	global $contigSeq;

	$unmappedContigs = '';
	$unmappedContigs_seq = '';

	$draftContig = array();
	$handle = fopen($draftFileName, 'r');
	while(($line = fgets($handle)) !== false){
		if($line[0] == '>'){
			array_push($draftContig, $line);
		}
	}
	fclose($handle);
	$n = 0;
	foreach($draftContig as $eachInDraft){
		$toBeTok = $eachInDraft;
		$contigNameInD = strtok(substr($toBeTok, 1), " \n");
		if($contigNameInD == '') continue;
		foreach($originMarkers as $eachInContig){
			$contigNameInC = strtok($eachInContig, " \n");
			if($contigNameInC == '') continue;
			if($contigNameInD == $contigNameInC){
				$unmapped = 0;
				break;
			}
			else{
				$unmapped = 1;
			}
		}
		if($unmapped == 1){
			$unmappedContigs .= "$contigNameInD $plus_o\n";
		}
	}
	file_put_contents($outputFile_draft, $unmappedContigs, FILE_APPEND);
}

function judge_dir($coords){
	$arr = file($coords, FILE_IGNORE_NEW_LINES);
	for($i = 0; $i < 5; $i++){
		array_shift($arr);
	}
	$element = preg_split('/ +/', trim($arr[0]));
	$num = count($element);
	$lines = array();
	foreach($arr as $k => $line){
		$element = preg_split('/ +/', trim($line));
		$lines[$k][7] = $element[7]; //que_len
		$lines[$k][$num-3] = $element[$num-3]; //ref_frm
		$lines[$k][$num-2] = $element[$num-2]; //que_frm

	}
	$line_num = count($lines);
	$forward_len = 0;
	$reverse_len = 0;
	for($i = 0; $i < $line_num; $i++){
		if($lines[$i][$num-3] * $lines[$i][$num-2] > 0){
			$forward_len += $lines[$i][7];
		}
		else{
			$reverse_len += $lines[$i][7];
		}
	}
	if($forward_len > $reverse_len){
		return 1;
	}
	else{
		return -1;
	}
}

function flip_scaffold($scaffold){
	if($scaffold[0][0] == '>'){
		$scaffold_name = $scaffold[0];
		array_shift($scaffold);
	}

	$scaffold = array_reverse($scaffold);
	foreach($scaffold as $k => $contig){
		list($contig_name, $ori) = preg_split('/ +/', $contig);
		$scaffold[$k] = $contig_name.' '.(1-intval($ori));
	}

	if(isset($scaffold_name))
		array_unshift($scaffold, $scaffold_name);

	return $scaffold;
}

function cmp($a, $b){
	if($a[0] == $b[0]){
		return ($a[1] < $b[1]) ? -1 : 1;
	}
	return ($a[0] < $b[0]) ? -1 : 1;
}

function cmp_arrSize($a, $b){
	$a_size = count($a);
	$b_size = count($b);
	if($a_size == $b_size){
		return 0;
	}
	return ($a_size < $b_size) ? 1 : -1;
}

function checkFASTA($file, $subject){
	global $meg, $meg_h;
	$handle = fopen($file, 'r');
	$contents = fread($handle, 1);
	if($contents != '>' && $contents != ';'){
		fclose($handle);
		$subject = $subject == 't' ? 'Target' : 'Reference';
		pErr("$subject is not a FASTA or multi-FASTA file");
	}
	fclose($handle);
}

function pErr($str){
	exit("ERROR: $str\n");
}

function p($str){
	echo "$str\n";
}

function getMicrotime(){
	list($usec, $sec) = explode(' ', microtime());
	return ((float)$usec + (float)$sec);
}
?>
