#! /usr/bin/python

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "gelpi"
__date__ = "$27-jun-2019 13:53:56$"

import re
import sys

#// Constants for filtering PDB coordinates
FIRSTHET = 1
GROUP = 2
NOHEAD = 3
HEAD = 4
ATOMSET = 5
NOALT = 6
NOWAT = 7
NOANISOU = 8
#// Groups
groupAtDefs = {
    'POLAR' :'N*,O*',
    'APOLAR' : 'C*,S*',
    'NOH' : 'C*,O*,S*,N*,P*',
    'BACK' : 'CA,O,N,NH,HA*',
    'NABACK': "C1',C2',C3',C4',C5',O1',O2',O3',O4',O5',P,OP1,OP2,H1',H2',H2'',H3',H4',H5',H5'',HO5'"
}
groupResDefs = {
    'POLAR' : 'ASN,ASP,GLN,GLU,SER,THR,TYR,HIS,LYS,ARG',
    'APOLAR' : 'ALA,VAL,LEU,ILE,MET,PRO,CYS,PHE,TRP',
    'NUC' : 'DA,DC,DG,DT,A,C,G,U',
    'PROT' : 'ASN,ASP,GLN,GLU,SER,THR,TYR,HIS,LYS,ARG,ALA,VAL,LEU,ILE,MET,PRO,CYS,PHE,TRP,GLY',
    'DNA' : 'DA,DC,DG,DT',
    'RNA' : 'A,C,G,U'
}

def filter(fmt, mode, pdb, filter='', invers=False): 
    if mode ==  FIRSTHET:
        ress = ""
        outpdb = ""
        for lin in pdb.split('\n'): 
            if re.match("/^HETATM/", lin) and re.match("/" + filter + "/", lin):
                if fmt == 'cif':
                    pdbData = re.split("/\s+/", lin)
                    resId = pdbData[6] + pdbData[7]
                elif fmt == 'pdb':
                    resId = lin[21:26]
                else:
                    print(fmt)
            if not ress:
                ress = resId
            if ress == resId:
                outpdb += lin+"\n"
    elif mode == HEAD:
        outpdb = headers(pdb, True)

    elif mode == NOHEAD:
        outpdb = headers(pdb, False)

    elif mode == NOANISOU:
        outpdb = ""
        filter = filter.upper()
        for lin in pdb.split('\n'):
            if not re.match('/^ANISOU/', lin):
                outpdb += lin + "\n"

    elif mode == GROUP: 
        outpdb = ""
        filter = filter.upper()
        negate = re.match('/^!/', filter)
        filter = filter.replace('!', '')
        filter = filter.replace('HETATM', 'HETATM,CONECT')
        regex = "/^(" + filter.replace(',', '|') + ")/"
        pdb = headers(pdb, False);
        for lin in pdb.split("\n"):
            if re.match('/^(MODEL|ENDMDL)/', lin):
                outpdb += lin + "\n"
                if not negate:
                    if re.match(regex, lin) :
                        outpdb += lin + "\n"
                else:
                    if not re.match(regex, lin):
                        outpdb += lin + "\n"

    elif mode == NOALT:
#            $outpdb = "";
#            if (!$filter) { // allows to select a specific alternative
#                $filter = "A";
#            }  
        if fmt == 'cif':
#            $filt = "\.$filter";

        elif fmt == 'pdb':
#                $filt = " $filter";
#            }
#            foreach (explode("\n", $pdb) as $lin) {
#                switch ($fmt) {
#                    case 'cif':
#                        $pdbData  = preg_split("/\s+/",$lin);
#                        $altId = $pdbData[4];
#                        break;
#                    case 'pdb':
#                        $altId = substr($lin, 16, 1);
#                        break;
#                }
#                if (!preg_match('/^(ATOM|HETATM)/', $lin) or preg_match('/[ ' . $filt . ']/', $altId)) {
#                    $outpdb .= "$lin\n";
#                }
#            }
#
#            break;
    elif mode == NOWAT:
#            $outpdb = "";
#            foreach (explode("\n", $pdb) as $lin) {
#                switch ($fmt) {
#                    case 'pdb':
#                        $resName = substr($lin, 17, 3);
#                        break;
#                    case 'cif':
#                        $pdbData  = preg_split("/\s+/",$lin);
#                        $resName = $pdbData[5];
#                        break;
#                }
#                
#                if (!preg_match('/(HOH|WAT)/', $resName)) {
#                    $outpdb .= "$lin\n";
#                }
#            }
#
#            break;
    elif mode == ATOMSET:
#//  [NomRes]NumRes:Chain.nomat/model Jmol like
#            /*
#             * ATOM     12  O   PRO A  47      47.248  16.227  63.502  1.00 61.58           O
#             * 012345678901234567890123456789012345678901234567890123456789012345678901234567
#             */
#            $outpdb = '';
#            $fd = parse_filter($filter);
#            $pdbnh = headers($pdb, False);
#            $model=0;
#            foreach (explode("\n", $pdbnh) as $lin) {
#                switch ($fmt) {
#                    case 'pdb':
#                        $resName  = substr($lin, 17, 3);
#                        $resId    = substr($lin, 22, 4);
#                        $chainId  = substr($lin, 21, 1);
#                        $atomName = substr($lin, 12, 4);
#                        break;
#                    case 'cif':
#                        $pdbData  = preg_split("/\s+/", $lin);
#                        $resName  = $pdbData[5];
#                        $resId    = $pdbData[8];
#                        $chainId  = $pdbData[6];
#                        $atomName = $pdbData[3];
#                        $model    = $pdbData[25];
#                        break;
#                }
#                $linok = false;
#                if (preg_match('/^MODEL *([0-9]*)/', $lin, $match)) {
#                    $model = $match[1];
#                }
#                if (preg_match('/^(MODEL|ENDMDL)/', $lin) and ( !$fd[5] or matchList($fd[5], $model))) {
#                    $linok = true;
#                } elseif (
#                        (!$fd[5] or ! $model or matchList($fd[5], $model)) and 
#                        (!$fd[1] or matchList($fd[1], $resName)) and 
#                        (!$fd[2] or matchList($fd[2], $resId))and 
#                        (!$fd[3] or matchList($fd[3], $chainId)) and 
#                        (!$fd[4] or matchList($fd[4], $atomName))
#                ) {
#                    $linok = true;
#                }
#                if ((!$not and $linok) or ( $not and ! $linok)) {
#                    $outpdb .= "$lin\n";
#                }
#            }
#            break;
    return outpdb
#}
#
#function parseAtGroup($groupAt) {
#    $atomset = [];
#    $groupAt = str_replace(' ', '', $groupAt);
#    foreach (explode(',', $groupAt) as $gr) {
#        $fixgr = preg_replace('/^!/', '', strtoupper($gr));
#        if (isset($GLOBALS['groupAtDefs'][$fixgr])) { //predefined groups
#            $atomset = array_merge($atomset, explode(',', $GLOBALS['groupAtDefs'][$fixgr]));
#        } else { // atoms
#            if (preg_match('/!/', $gr)) {
#                //TODO
#            } else {
#                $atomset[] = $gr;
#            }
#        }
#    }
#    return ':.' . join(",", $atomset);
#}
#
#function parse_filter($a) {
#    $d = [];
#    $a = str_replace(' ', '', $a);
#    #$a = str_replace('*', '', $a);
#// afegir camps buits
#    if (!preg_match('/\[/', $a)) {
#        $a = '[]' . $a;
#    }
#    if (!preg_match('/:/', $a)) {
#        $a = preg_replace('/(\][0-9\-\,]*)/', "$1:", $a);
#    }
#    if (!preg_match('/\./', $a)) {
#        $a = preg_replace('/:([^\/]*)/', ":$1.", $a);
#    }
#//
#    $d = preg_split('/[\[\]:\.\/]/', $a);
#// Mantenim * nomes a atom
#    $d[2] = str_replace('*', '', $d[2]);
#    $d[3] = str_replace('*', '', $d[3]);
#    if (count($d)< 6) { // Avoid warnings for incomplete filters
#        $d[5]='';
#    }
#    return $d;
#}
#
#function headers($pdb, $rethead) {
#    $outpdb = "";
#    $head = True;
#    foreach (explode("\n", $pdb) as $lin) {
#        if (preg_match("/^(ATOM|HETATM|MODEL)/", $lin)) {
#            $head = False;
#        }
#        if (preg_match("/^#/", $lin)) {
#            $head = True;
#        }
#        if ($head == $rethead) {
#            $outpdb .= "$lin\n";
#        }
#    }
#    return $outpdb;
#}
#
#function matchList($l, $a) {
#    $m = false;
#    foreach (explode(',', $l) as $t) {
#        if (preg_match('/-/', $t)) {
#            list ($min, $max) = explode('-', $t);
#            $m = ($m or ( (trim($a) >= $min) and ( trim($a) <= $max)));
#        } else {
#            $m = ($m or matchStr(trim($t), trim($a)));
#        }
#    }
#    return $m;
#}
#
#function matchStr($a, $b) {
#    return (fnmatch($a, $b) or fnmatch($b, $a));
#}
#
#function calcHBonds($pdb, $inter = false, $type = false) {
#    $pdbnh = headers($pdb, False);
#    $outpdb = '';
#    header('Content-type:text/plain');
#    $i = 0;
#    foreach (explode("\n", $pdbnh) as $lin) {
#        $nat[$i] = substr($lin, 5, 6) + 0;
#        $atid[$i] = trim(substr($lin, 12, 5));
#        $resid[$i] = trim(substr($lin, 17, 4));
#        $chid[$i] = substr($lin, 21, 2);
#        $resn[$i] = substr($lin, 23, 4);
#        $residstr[$i] = $resid[$i] . ' ' . $chid[$i] . trim($resn[$i]);
#        $atx[$i] = substr($lin, 30, 8) + 0;
#        $aty[$i] = substr($lin, 38, 8) + 0;
#        $atz[$i] = substr($lin, 46, 8) + 0;
#        if (preg_match('/^ATOM|HETA/', $lin)) {
#            $i++;
#        }
#    }
#    $totat = $i - 1;
#    for ($i = 0; $i < $totat - 1; $i++) {
#        $hbs[$i] = [];
#        for ($j = $i + 1; $j < $totat; $j++) {
#            if (($residstr[$i] != $residstr[$j]) and ( !$inter or ( $chid[$i] != $chid[$j]))) {
#                $d = calcDist2($atx[$i], $aty[$i], $atz[$i], $atx[$j], $aty[$j], $atz[$j]);
#                if ($d <= $GLOBALS['HBDIST2']) {
#                    $hbs[$i][] = [$j, $d];
#                }
#            }
#        }
#        foreach (array_values($hbs[$i]) as $hb) {
#            $tipHB = $GLOBALS['HBondDef'][$resid[$i] . ' ' . $atid[$i] . ' ' . $resid[$hb[0]] . ' ' . $atid[$hb[0]]] . $GLOBALS['HBondDef'][$resid[$hb[0]] . ' ' . $atid[$hb[0]] . ' ' . $resid[$i] . ' ' . $atid[$i]];
#            if (!$type or ( $type == $tipHB)) {
#                $outpdb .= "$residstr[$i] $atid[$i] - " . $residstr[$hb[0]] . " " . $atid[$hb[0]] . " : " . sprintf("%5.2f", sqrt($hb[1]));
#            }
#            if ($tipHB) {
#                $outpdb .= "($tipHB)";
#            }
#            $outpdb .= "\n";
#        }
#    }
#    return $outpdb;
#}
#
#function calcDist2($x1, $y1, $z1, $x2, $y2, $z2) {
#    return ($x1 - $x2) * ($x1 - $x2) + ($y1 - $y2) * ($y1 - $y2) + ($z1 - $z2) * ($z1 - $z2);
#}
#
#function searchPDB($params) {
#    $cond = [['supersby' => ['$exists' => False]]];
#    if (isset($params['sequence'])) {
#        if ($params['seqType'] == "exact") {
#            $seqc = $sequencesCol->find([
#                'sequence' => $params['sequence'],
#                'origin' => 'pdb',
#                'type' => $params['molTy']], ['_id' => 1]);
#            $seqc->timeout(-1);
#        } else {
#            $seqQuery = str_replace('X', '.', strtoupper($params['sequence']));
#            $seqQuery = str_replace('-', '', $seqQuery);
#            $seqQuery = str_replace('(', '{', $seqQuery);
#            $seqQuery = str_replace(')', '}', $seqQuery);
#            $seqregex = new MongoRegex("/" . $seqQuery . "/");
#            $conds = ['sequence' => $seqregex, 'origin' => 'pdb'];
#            if ($params['molTy'] != 'Any') {
#                $conds['type'] = $params['molTy'];
#            }
#//        print "<pre>";
#//        print json_encode($conds);
#//        print "</pre>";
#            $seqc = $GLOBALS['sequencesCol']->find($conds, ['_id' => 1]);
#            $seqc->timeout(-1);
#        }
#        foreach (array_keys(iterator_to_array($seqc)) as $sqid) {
#            $ids[substr($sqid, 0, 4)] = 1;
#        }
#        if (!count($ids)) {
#            return [];
#        }
#        $cond[] = ['_id' => ['$in' => array_keys($ids)]];
#    }
#    if (isset($params['resmin']) or isset($params['resmax'])) {
#        if (isset($params['resmin'])) {
#            $cond[] = ['resol' => ['$gt' => $params['resmin']]];
#        }
#        if (isset($params['resmax'])) {
#            $cond[] = ['resol' => ['$lt' => $params['resmax']]];
#        }
#    }
#    if (!isset($params['anyComp'])) {
#        $cl2 = [];
#        foreach (array_keys($params['compType']) as $ct) {
#            $cl2[] = ['compType' => $ct];
#        }
#        if (count($cl2) > 1) {
#            $cond[] = ['$or' => $cl2];
#        } else {
#            $cond[] = ['compType' => $cl2[0]['compType']];
#        }
#    }
#    if (!isset($params['anyExpType'])) {
#        $cl2 = [];
#        foreach (array_keys($params['expType']) as $ct) {
#            $cl2[] = ['expType' => $ct];
#        }
#        if (count($cl2) > 1) {
#            $cond[] = ['$or' => $cl2];
#        } else {
#            $cond[] = ['expClasse' => $cl2[0]['expType']];
#        }
#    }
#    if (isset($params['query'])) {
#// Text index pendents versio Mongodb
#        foreach (explode(' ', $params['query']) as $wd) {
#            $cl2 = [];
#            foreach (array_keys($params['queryOn']) as $fld) {
#                $rex = new MongoDB\BSON\Regex($wd, "i");
#                $cl2[] = [$fld => $rex];
#            }
#            if (count($cl2) > 1) {
#                $cond[] = ['$or' => $cl2];
#            } else {
#                $cond[] = $cl2[0];
#            }
#        }
#    }
#    if (count($cond)) {
#        $fcond = ['$and' => $cond];
#    } else {
#        $fcond = [];
#    }
#//print "<pre>";
#//print json_encode($fcond);
#//print "</pre>";
#// Fasta output (experimental)
#    if (isset($params['fmt']) and ($params['fmt'] == 'fasta')) {
#        $filter = ['_id' => 1, 'chain' => 1];
#        $chainIds = [];
#        foreach (iterator_to_array($GLOBALS['PDB_EntryCol']->find($fcond, $filter)) as $ch) {
#            $chainIds = array_merge($chainIds, $ch['chain']);
#        }
#        $results = $GLOBALS['sequencesCol']->find(['_id' => ['$in' => $chainIds]]);
#    } else {
#        $filter = [];
#        foreach (explode(',', str_replace(' ', '', $params['fields'])) as $f) {
#            $filter[$f] = 1;
#        }
#        if (isset($params['sort'])) {
#            switch ($params['sort']) {
#                case 'header': $sortA = ['header' => 1, '_id' => 1];
#                    break;
#                case 'compType': $sortA = ['compType' => 1, '_id' => 1];
#                    break;
#                case 'expType': $sortA = ['expType' => 1, '_id' => 1];
#                    break;
#                case 'resol': $sortA = ['resol' => 1, '_id' => 1];
#                    break;
#                case '_id':
#                default:
#                    $sortA = ['_id' => 1];
#            }
#        } else {
#            $sortA =[];
#        }
#        $resultCount = $GLOBALS['PDB_EntryCol']->count($fcond);
#        $results = $GLOBALS['PDB_EntryCol']->find($fcond, ['projection' => $filter,'sort'=> $sortA]);
#    }
#    return ["search"=>$results,"count"=>$resultCount];
#}
#
#function getPDBwithMon($idMon) {
#    $data = $GLOBALS['PDB_EntryCol']->findOne(['hetAtms' => $idMon], ['_id' => 1]);
#    return $data['_id'];
#}
#
#function getMonomerData($idMon, $pdbs = false) {
#    $data = $GLOBALS['monomersCol']->findOne(['_id' => $idMon]);
#    if ($pdbs) {
#        $pdbs = iterator_to_array(
#                $GLOBALS['PDB_EntryCol']->find(
#                        ['hetAtms' => $idMon],
#                        [
#                            'projection' => ['_id' => 1],
#                            'sort'=> ['_id' => 1]
#                        ]));
#        $data['pdbs']=[];
#        foreach ($pdbs as $p) {
#            $data['pdbs'][] = $p['_id'];
#        }
#    }
#    return $data;
#}
#
#/* mmCIF support
# * parsemmCif generate a PHP array structure, that can be shown as json or xml 
# */
#
#function parsemmCif($cifData, $params) {
#    set_time_limit(0);
#    $groups = explode("#", $cifData);
#    $data = [];
#    while ($gr = array_shift($groups)) {
#        $data = array_merge($data, _parsemmCifGroup($gr, $params));
#    }
#    return $data;
#}
#
#function _parsemmCifGroup($gr, $params) {
#    $lines = explode("\n", $gr);
#    array_shift($lines);
#    if ($lines[0] == "loop_") {
#        array_shift($lines);
#        return _parsemmCifGroupLoop($lines, $params);
#    } else {
#        return _parsemmCifGroupNoLoop($lines, $params);
#    }
#}
#
#function _parsemmCifGroupNoLoop($lines, $params) {
#    $data = [];
#    $nl = 0;
#    $v = '';
#    while ($nl < count($lines)) {
#        $l = _protectQuotedBlock($lines[$nl]);
#        if (preg_match("/^_/", $l)) {
#            $data = _addValue($data, $k, _unprotectQuotedBlock($v), $params);
#            list ($k, $v) = preg_split("/\s+/", $l);
#        } else {
#            $v .= preg_replace("/^;/", "", $l);
#        }
#        $nl++;
#    }
#    $data = _addValue($data, $k, _unprotectQuotedBlock($v), $params);
#    return $data;
#}
#
#function _parsemmCifGroupLoop($lines, $params) {  
#    $data = [];
#    $fields = [];
#    foreach ($lines as $l) {
#        if (preg_match("/^_/", $l)) {
#            list ($grup, $f) = explode(".", preg_replace("/^_/", "", $l));
#            $fields[] = $f;
#        }
#    }
#    if (!preg_match("/s$/", $grup) and ( $params['fmt'] == 'xml')) { // management lists in xml through plural/singular tags
#        $grup .= "s";
#        $grup = preg_replace("/ys$/", "ies", $grup); // proper plural when ..y :-)
#    }
#    if (preg_match ("/atom_site/",$grup)) { 
#        switch ($params['atoms']) {
#            case "full":
#                break;
#            case "none":
#            case "no":
#                return [];
#                break;
#            case "block":
#            default:
#                $data[$grup]['fieldList'] = join (" ", $fields);
#                $fields = ['atom_records'];
#                $nl =0;
#                while ($nl < count($lines)) {
#                    if (!preg_match("/^_/", $lines[$nl])) {
#                        $data[$grup]['atom_records'][] = $lines[$nl];
#                    }
#                    $nl++;
#                }
#                return $data;
#        }
#    }
#    $nl = 0;
#    $dataFields = [];
#    while ($nl < count($lines)) {
#        if (!preg_match("/^_/", $lines[$nl])) {
#            if (preg_match("/^;/", $lines[$nl])) { // read block between ";" and add as a quoted single line 
#                $concatBl = $lines[$nl];
#                $nl++;
#                while ($nl < count($lines) and ! preg_match("/^;/", $lines[$nl])) {
#                    $concatBl .= $lines[$nl];
#                    $nl++;
#                }
#                $concatBl .= $lines[$nl];
#                $lines[$nl] = str_replace(";", "'", $concatBl);
#                $concatBl = '';
#            }
#            $l = _protectQuotedBlock($lines[$nl]);
#            $newData = preg_split("/\s+/", $l);
#            $dataFields = array_merge($dataFields, $newData);
#
#            }
#        $nl++;
#    }
#    // Distribute $dataFields among n records of m fields each, skip empties
#    $nf = 0;
#    while ($nf < count($dataFields)) {
#        $dataList = [];
#        for ($i = 0; $i < count($fields); $i++) {
#            if ($dataFields[$nf] === '') {
#                $nf++;
#            }
#            $dataList = _addValue($dataList, $fields[$i], _unprotectQuotedBlock($dataFields[$nf]), $params);
#            $nf++;
#        }
#        if ($dataList) {
#            $data[$grup][] = $dataList;
#        }
#    }
#    return $data;
#}
#
#function _addValue($data, $k, $v, $params) {
#    $k = preg_replace("/^_/", "", preg_replace("/\s+$/", "", $k));
#    if (isset($params['noEmpty'])) {
#        $v = preg_replace("/^\?$/", "", $v);
#    }
#    if ($k and ( $v !== '')) {
#        if (!preg_match('/\./', $k)) {
#            if (($params['fmt'] == "xml")) {
#                $k = preg_replace("/s$/", "s_", $k); // avoid open lists in xml due to final "s"
#            }
#            $data[$k] = $v;
#        } else {
#            $kk = explode(".", $k);
#            if (($params['fmt'] == "xml")) {
#                $kk[0] = preg_replace("/s$/", "s_", $kk[0]);
#                $kk[1] = preg_replace("/s$/", "s_", $kk[1]);
#            }
#            $data[$kk[0]][$kk[1]] = $v;
#        }
#    }
#    return $data;
#}
#
#function _protectQuotedBlock($l) { // process 'ccc, cc' en 1 block
#    return preg_replace_callback(
#            "/'([^']*)'/", 
#            function ($a) {
#                return preg_replace("/\s/", '#', $a[0]);
#            }
#            , $l);
#}
#
#function _unprotectQuotedBlock($l) {
#    return str_replace("#", " ", str_replace("'", "", $l));
#}
#
#function getPDBPict ($params) {
#// TODO $bunit
#    if (isset($params['chain'])) {
#       $fn = strtoupper($params['id'])."_".$params['chain'];
#    } else {
#       $fn = strtoupper($params['id']);
#    }
#    $ff = getGSFile($GLOBALS['PDB_PictsGrid'], "$fn.png");
#    return $ff;
#}


if __name__ == "__main__":
    pdb_input_path = sys.argv[0]
    pdb_output_path = sys.argv[1]
    operation = sys.argv[2]
    
