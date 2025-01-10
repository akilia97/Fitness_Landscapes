#!/bin/bash
# Purpose: To Visualize influenza virus evolutionary trajectories in multidimensional landscapes
# Author: Akilia Mathie
# Last Modified: 2024-07-22

### DATABASE VARIABLES ###
hDB=protein_modeling #this is the schema that you want the sql query to run or where you want to pull/push info from
hPATH=/warehouse/tablespace/external/hive/${hDB}.db #files live on hive

### DEFINE BASE PATH ###
bpath=
if [ "$bpath" == "" ]; then
    bpath=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi

### LOAD DIRECTORY PATHS ###
binp=$bpath/bin
libp=$bpath/lib
sqlp=$bpath/sql
tblp=$bpath/tbl

### DEFINE FUNCTIONS ###
PROGRAM=FitnessLandscape

function time_stamp() {
	local t=$(date +"%Y-%m-%d %k:%M:%S")
	echo -e "[$t]\t$PROGRAM ::: $1"
}

function die() {
	local t=$(date +"%Y-%m-%d %k:%M:%S")
        echo -e "[$t]\t$PROGRAM ERROR :: $1" 1>&2
        exit 1
}

function print_help {
    cat << EOL

    Valid input format: ./$PROGRAM <FUNCTION>
       
    Functions:  autoML [INPUT] [OUTPUT] [DELIMITER] [COMPOSITE|MIXED] [ALL|ANTIGENIC|PROTEOMIC]
                distance
                do
                impute
                merge
                predictML [INPUT] [OUTPUT] [DELIMITER] [ITERATIVE] [ALL|ANTIGENIC|MIXED|PROTEOMIC] [REFERENCE INPUT]
                update

EOL
    exit 0
}

#waiting on potentially a 3rd sql query for structural distance

## RUN PIPELINE FROM START TO FINISH ##
if [ "$1" == "do" ]; then
    time_stamp "BEGIN"

    # Step 1: Update local files with antigenic, gentic, and proteomic trends hosted on NDW Cloudera Data Platform
    # if [[ -f "tbl/structural_distance_calpha.txt" ]]; then
    #    vhs=$( awk -F "\t" 'NR != 1 { print $4 ":" }' "tbl/structural_distance_calpha.txt" | sort -u )
    # else 
    #    vhs=()
    # fi
    
    # vhsd="$( echo ${vhs} )"
    # vhsf=$( echo $vhsd | sed 's/ //g' | sed 's/.$//' )
    
    sql_queries=( indv_antigenic_summaries antigenic_summaries genetic_summaries rosetta_summaries )
    for sql in "${sql_queries[@]}"; do
        # cp $sqlp/templates/${sql}_summaries_template.sql $sqlp/${sql}_summaries.sql
        # sed -i "s|INSERTVARIANTHASHLIST|${vhsf}|g" $sqlp/${sql}_summaries.sql
        # sed -i "s|\:|\'\,\n\'|g" $sqlp/${sql}_summaries.sql

        # Step 1A: Execute 'himpala' command to download antigenic/genetic trends from NDW CDP tables
        himpala "$sqlp/${sql}.sql" -p > "$tblp/${sql}.txt" || die "Extract of '${sql}.sql' failed. Exiting."
        time_stamp "Extracted ${sql} trends"
    done

    # Step 1B: Execute context-specific structural distance (SD) SQL queries via the 'himpala' command
    # REFRESH local context-specific SD tab-delimited TXT files
    contexts=( calpha carbon heavy )
    for context in "${contexts[@]}"; do
        cp $sqlp/templates/structural_distance_template.sql $sqlp/${context}_structural_distance.sql
        sed -i "s|ATOMICCONTEXT|${context}|g" $sqlp/${context}_structural_distance.sql

        himpala "$sqlp/${context}_structural_distance.sql" -p > "$tblp/structural_distance_${context}.txt" || die "Extract of '${context}_structural_distance.sql' failed. Exiting."
        time_stamp "Extracted ${context} SD metrics"
    done

    # Step 1C: Refresh local lot-specific antisera data dictionary
    himpala "$sqlp/reportable_antiserum.sql" -p > "$tblp/seasonal_antisera_references.txt" || die "Extract of 'reportable_antiser
um.sql' failed. Exiting."

    # Step 2: Merge all local genetic, proteomic and antigenic data
    time_stamp "Merged all local data files."
    python $binp/data_merge.py || die "Data merge failed"
    python $binp/data_merge.py -i True || die "Data merge (individual) failed"
    
    # Step 3: Impute (or replace) "missing" data (i.e., NaNs) in provided data matrices
    subtypes=( B H1 H3 )
    for subtype in "${subtypes[@]}"; do
        time_stamp "Imputed ${subtype} data/metadata."
        python "$binp/data_imputation.py" -i "$tblp/${subtype}_consolidated.txt" -d tab -o "$tblp/${subtype}_iconsolidated.txt" || die "Data imputation failed for ${subtype}"
        python "$binp/data_imputation.py" -i "$tblp/indv_${subtype}_consolidated.txt" -d tab -o "$tblp/indv_${subtype}_iconsolidated.txt" || die "Data imputation failed for ${subtype} (individual)"
    done

    # Step 4: Implementate ensemble machine-learning (ML) methods to compute "predicted" landscape for applied XYZ factor coordinates
    # NOTE: Need to develop an additional BASH script to submit to SCICOMP HPC
    subtypes=( B H1 H3 )
    for subtype in "${subtypes[@]}"; do
        python "$binp/predictML.py -i $tblp/${subtype}_iconsolidated.txt -o $tblp/${subtype}_plandscape.txt -d tab -z antigenic -b $tblp/seasonal_antisera_references.txt" || die "predictML failed for $subtype"
    done

    time_stamp "END"

## IDENTIFY OPTIMAL ML MODELS FOR USER-DEFINED DATASETS ##
elif [ "$1" == "autoML" ]; then
    time_stamp "BEGIN"
    
    python "$binp/autoML.py" -i $1 -o $2 -d $3 -x $4 -y $5 || die "autoML failed"

    time_stamp "END"

## COMPUTE STRUCTURAL DISTANCE METRIC ##
elif [ "$1" == "distance" ]; then
    time_stamp "BEGIN"

    contexts=( calpha carbon heavy )
    for context in "${contexts[@]}"; do
        time_stamp "Updated ${context} atom structural distance metric(s)."

        if [[ -f "$tblp/structural_distance_${context}.txt" ]]; then 
            vhs=$( awk -F "\t" '{ print $4 }' "$tblp/structural_distance_${context}.txt" | sort -u )
        else 
            vhs=()
        fi
			
        for i in data/alphafold2/reference/*.pdb; do
            j=$(echo ${i} | awk -F "/" '{ print $4 }')
            
            if ! [[ ${vhs[*]} =~ "$j" ]]; then
                python $binp/structural_distance.py --input data/alphafold2/reference/${j} --output $tblp/structural_distance_${context}.txt --context ${context}
            fi
        done
    done        

    time_stamp "END"

## IMPUTE MISSING DATA IN LOCAL FILES ##
elif [ "$1" == "impute" ]; then
    time_stamp "BEGIN"

    subtypes=( B H1 H3 )
    for subtype in "${subtypes[@]}"; do
        time_stamp "Imputed ${subtype} data/metadata."
        python "$binp/data_imputation.py" -i "$tblp/${subtype}_consolidated.txt" -d tab -o "$tblp/${subtype}_iconsolidated.txt" || die "Data imputation failed for ${subtype}"
        python "$binp/data_imputation.py" -i "$tblp/indv_${subtype}_consolidated.txt" -d tab -o "$tblp/indv_${subtype}_iconsolidated.txt" || die "Data imputation failed for ${subtype} (individual)"
    done
    
    time_stamp "END"

## MERGE LOCAL DATA FILES ##
elif [ "$1" == "merge" ]; then
    time_stamp "BEGIN"
    
    python $binp/data_merge.py || die "Data merge failed"
    python $binp/data_merge.py -i true || die "Data merge (individual) failed"

    time_stamp "END"

## MODEL PREDICTED XYZ-COORDINATE LANDSCAPES ##
elif [ "$1" == "predictML" ]; then
    time_stamp "BEGIN"

    python "$binp/predictML.py" -i $2 -o $3 -d $4 -x $5 -z $6 -b $7 || die "predictML execution attempt failed"
    
    time_stamp "END"

## REFRESH LOCAL DATA FILES ##
elif [ "$1" == "update" ]; then
    time_stamp "BEGIN"

    sql_queries=( indv_antigenic_summaries antigenic_summaries genetic_summaries rosetta_summaries )
    for sql in "${sql_queries[@]}"; do
        himpala "$sqlp/${sql}.sql" -p > "$tblp/${sql}.txt" || die "Extract of '${sql}.sql' failed. Exiting."
        time_stamp "Extracted ${sql} trends"
    done

    contexts=( calpha carbon heavy )
    for context in "${contexts[@]}"; do
        cp $sqlp/templates/structural_distance_template.sql $sqlp/${context}_structural_distance.sql
        sed -i "s|ATOMICCONTEXT|${context}|g" $sqlp/${context}_structural_distance.sql

        himpala "$sqlp/${context}_structural_distance.sql" -p > "$tblp/structural_distance_${context}.txt" || die "Extract of '${context}_structural_distance.sql' failed. Exiting."
        time_stamp "Extracted ${context} SD metrics"
    done

    himpala "$sqlp/reportable_antiserum.sql" -p > "$tblp/seasonal_antisera_references.txt" || die "Extract of 'reportable_antiser
um.sql' failed. Exiting."

    time_stamp "END"

else
    print_help
    time_stamp "Review input commands. Entry invalid."
fi
