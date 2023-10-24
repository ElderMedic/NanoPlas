#!/bin/sh
usage(){
	echo "[usage] ./$0 [input.vcf] [output.csv]"
}

if [ $# -ne 2 ];then
	usage
	exit 1
fi

if [ $# -eq 2 ];then
	input_vcf=$1
	output=$2
	if [ ! -f $input_vcf ] ;then
		echo "$1 is non exist"
		exit 1
	fi

	grep -Ev "#" $input_vcf | awk -F "\t" 'BEGIN{print "Contig,Position,Ref,Alt,Type,Ref_COUNT,Alt_COUNT,TOTAL_COUNT"}
	{
	    split($8, dp, "DP4=");
	    split(dp[2], dp, ";");
	    split(dp[1], eachdp, ",");
	    split($5, gt, ",");
         tt = "single mismatch";

    if (length($4) > length($5)) {
        split($8, type, ";");
        tt = type[1] "[deletion]";
    } else if (length($4) < length($5)) {
        split($8, type, ";");
        tt = type[1] "[insertion]";
    } else if ($5 == ".") {
        tt = ".";
    }
	    #if(length($4)>length($5)){split($8,type,";");tt=type[1]"[deletion]"};
	    #if(length($4)<length($5)){split($8,type,";");tt=type[1]"[insertion]"};
	    #if(length($4)==length($5)){tt="single mismatch"};
	    printf $1","$2","$4","$5","tt",";
	    printf eachdp[1]+eachdp[2]",";
	    for(y=1; y<=length(gt); y++){
	        if(gt[y]!="<*>"){
	            altdp=gt[y];
                         split(altdp, alt_counts, ":");
                         printf alt_counts[2] + eachdp[3] + eachdp[4] ",";
	        }
	    }
	print eachdp[1]+eachdp[2]+eachdp[3]+eachdp[4];
	}' > $output
fi
