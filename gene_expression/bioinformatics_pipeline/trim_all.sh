#BSUB -u and128@miami.edu

and="/scratch/projects/and_transcriptomics"
cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.gz))

for sample in ${data[@]} ;
do \
echo '#!/bin/bash' > "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -q general' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -J '"${sample}"_trim_all'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -o '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"$sample"_trim_all%J.out'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -e '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"$sample"_trim_all%J.err'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -n 8' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '#BSUB -N' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'sample_name=$(basename "$sample" .fastq.gz)' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'if [[ "$sample" == *.gz ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '  gunzip -c "$sample" > "$sample_name.temp.fastq"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '  sample="$sample_name.temp.fastq"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'lead="${2:-[ATGC]?[ATGC][AC][AT][AT][AC][AT][ACT]GGG+|[ATGC]?[ATGC][AC][AT]GGG+|[ATGC]?[ATGC]TGC[AC][AT]GGG+|[ATGC]?[ATGC]GC[AT]TC[ACT][AC][AT]GGG+}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'trim=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'name=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'name2=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'seq=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'qua=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'declare -A seen' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'll=3' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'nohead=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'dups=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'ntag=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'tot=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'goods=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'goods=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'while IFS= read -r line; do' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    if [ "$ll" -eq 3 ] && [[ "$line" =~ ^@(.+)$ ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        tot=$((tot + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        name2="${BASH_REMATCH[1]}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        if [[ "$seq" =~ ^($lead)(.+) ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            start="${BASH_REMATCH[2]:0:20}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            idtag="${BASH_REMATCH[1]}$start"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            if [ -z "${seen[$idtag]}" ] && [[ "$idtag" != *N* ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                seen["$idtag"]=1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                trim="${#BASH_REMATCH[1]}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                echo "$name"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                echo "${BASH_REMATCH[2]}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                echo "+"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                echo "${qua:trim}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                goods=$((goods + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            elif [ -n "${seen[$idtag]}" ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                dups=$((dups + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            else' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '                ntag=$((ntag + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        elif [ "$tot" -gt 1 ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '            nohead=$((nohead + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        seq=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        ll=0' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        qua=""' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        name="$name2"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    elif [ "$ll" -eq 0 ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        seq="$line"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        ll=1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    elif [ "$ll" -eq 2 ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        qua="$line"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        ll=3' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    else' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        ll=2' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'done < "$sample".temp.fastq' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'if [[ "$seq" =~ ^($lead)(.+) ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    start="${BASH_REMATCH[2]:0:20}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    idtag="${BASH_REMATCH[1]}$start"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    if [ -z "${seen[$idtag]}" ] && [[ "$idtag" != *N* ]]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        seen["$idtag"]=1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        trim="${#BASH_REMATCH[1]}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        echo "$name"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        echo "${BASH_REMATCH[2]}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        echo "+"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        echo "${qua:trim}"' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        goods=$((goods + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    elif [ -n "${seen[$idtag]}" ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        dups=$((dups + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    else' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '        ntag=$((ntag + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'elif [ "$tot" -gt 1 ]; then' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '    nohead=$((nohead + 1))' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'fi' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'done > "$sample".trim.fastq' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'echo "$sample	total:$tot	goods:$goods	dups:$dups	noheader:$nohead	N.in.header:$ntag" > "$sample"_trimoutput.txt' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '${and}/programs/TrimGalore-0.6.10/trim_galore $sample.trim.fastq \' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '--adapter "AGATCGG" \' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '--adapter "AAAAAAAA" \' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '--quality 15 \' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '--length 25 \' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '--output_dir ${and}/Ch2_temperaturevariability2023/2_trimmed_reads/take_3 ;' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'done' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'cd ${and}/Ch2_temperaturevariability2023/2_trimmed_reads/take_3' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'for f in *.fq;' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'do' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'mv -v -- "$f" "${f%.fq}.fastq";' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'done' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo '' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'module load fastqc/0.10.1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'fastqc *.fastq' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
echo 'multiqc *' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job
bsub < "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/take_3/"${sample}"_trim_all.job ; \
done
