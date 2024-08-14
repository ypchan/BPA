#!/usr/bin/env bash

function help() {
	echo ""
	echo "convert_mrbayes_sampletrees_for_logcombiner.sh -- convert mrbayes sample trees file(*.run#.t) for the beast tool logcombiner"
	echo ""
	echo "DATE: 2021-03-06"
	echo "BUGS: Any bugs should be reported to chenyanpeng1992@outlook.com"
	echo "Example:"
	echo "    convert_mrbayes_sampletrees_style_for_logcombiner.sh mrbayes.run1.t mrbayes.run1.trees"
	echo "    convert_mrbayes_sampletrees_style_for_logcombiner.sh mrbayes.run2.t mrbayes.run2.trees"
	echo "    logcombiner -log mrbayes.run1.trees -log mrbayes.run2.trees -o mrbayes.trees"
	echo "    treeannotator -target mltree.newick mrbayes.trees mltree.bayesioninfo.newick"
	echo ""
	exit 1
}

function getHead() {
	echo $1 $2
    head -n 1 $1 > $2
    sed -i '1G' $2
}

function addTaxaBlock() {
    taxonNum=$(grep -e '^      ' $1 | wc -l)
    echo "Begin taxa;" >> $2
    echo "    Dimensions ntax=${taxonNum};" >> $2
    echo "        Taxlabels" >> $2
	grep '^      ' $1 | awk '{ print "            "$2 }' | sed 's/,//;s/;//' >> $2
    echo '            ;' >> $2
    echo 'End;' >> $2
}

function modifyTreeBlock() {
    echo "Begin trees;" >> $2
    echo "    Translate" >> $2
	grep '^      ' $1 | sed 's/^      /        /;s/;//' >> $2
	echo ';' >> $2
	grep '^   tree gen' $1 | sed 's/^   //;s/gen\./STATE_/;s/\[&U\] // ' >> $2
	echo 'Ends;' >> $2
}

# Main
if [[ $# != 2 ]]; then
    help
fi

getHead $1 $2
addTaxaBlock $1 $2
modifyTreeBlock $1 $2
