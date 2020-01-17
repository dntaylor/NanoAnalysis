nCores=20


function waitForJobs {
    joblist=($(jobs -p))
    while (( ${#joblist[@]} >= $nCores ))
    do
        sleep 1
        joblist=($(jobs -p))
    done
}

function finishJobs {
    joblist=($(jobs -p))
    while (( ${#joblist[@]} >= 0 ))
    do
        sleep 1
        joblist=($(jobs -p))
    done
}

for proc in processors/*.py; do
    procpy=`basename $proc`
    baseProcessor="${procpy%.*}"
    for year in 2016 2017 2018; do
        processor=processors/${baseProcessor}_${year}.coffea
        for fileset in filesets/$year/*.json; do
            dataset=$(basename "$fileset" .json)
            if [[ "$dataset" == "all" ]]; then
                continue
            fi
            if [[ "$dataset" == "data" ]]; then
                continue
            fi
            if [[ "$dataset" == "mc" ]]; then
                continue
            fi
            hists=hists/${baseProcessor}/${year}/${dataset}.coffea
            if [ $hists -nt $processor ] && [ $hists -nt $fileset ]; then
                #echo "Skipping" $baseProcessor $year $dataset
                continue
            fi
            nfiles=`jq -r ".[\"${dataset}\"]" $fileset | jq length`
            njobs=$((nfiles/2+1))
            waitForJobs
            #echo "Launching" $baseProcessor $year $dataset
            ./run_processor.py -j $njobs -r 1 --parsl --condor $baseProcessor $year $fileset &
            # sleep for a while so that the new parsl run can start
            sleep 20
        done
    done
done
