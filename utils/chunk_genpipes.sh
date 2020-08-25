#!/usr/bin/env bash

# FUNCTIONS

usage (){

echo
echo "usage: $0 <GENPIPES SCRIPT> <OUTPUT FOLDER>
  Chunk Genpipes submiting sript so there is njobs in them"
echo
echo "   <GENPIPES SCRIPT>       A Genpipes output script"
echo "   <OUTPUT FOLDER>         Folder where to store chunks"
echo "   -n <chunk size>         Maximum number of job in chunk,"
echo "                             default=20"

}

load_previous_submit_id (){
echo load_previous_submit_id on $1
 current_file=$1
cat << 'EOF' >> ${current_file}
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
for file in $(ls ${SCRIPTPATH}/*out | sort -n ); do
source $file
done
EOF
}

start_new_chunk () {
echo start_new_chunk on $1
  cat << EOF > $1
#!/usr/bin/env bash
source $2
STEP=$3
mkdir -p \$JOB_OUTPUT_DIR/\$STEP
EOF
}

# SCRIPT


slurm_user=$USER

max_chunk=20
while getopts "hn:u:" opt; do
  case $opt in
    u)
      slurm_user=${OPTARG}
    ;;
    n)
      max_chunk=${OPTARG}
    ;;
    h)
      usage
      exit 0
    ;;
    \?)
      usage
      exit 1
    ;;
  esac
done

shift $((OPTIND-1))


if [ $# -lt 2 ]; then
  usage
  exit 1
fi
genpipes_in=$1
out_dir=$2



mkdir -p ${out_dir}
TIMESTAMP=`date +%FT%H.%M.%S`
header=${out_dir}/header.sh
echo '# header for all chunks' > ${header}
while read -r line ; do
    if [[ $line =~ ^STEP=.*$ ]]; then
      STEP=${line#STEP=}
      break
    elif [[ "$line" =~ ^TIMESTAMP=.*$ ]]; then
      line="TIMESTAMP=${TIMESTAMP}"
    fi
    echo "$line" >> ${header}
done < ${genpipes_in}
echo "CHUNK_SIZE=${max_chunk}" >> ${header}

chunk=1
out_file=/dev/null
n_job=0
# create chunks
while read -r line ; do
    echo "$line" >> $out_file

    if [ "$line" == 'cd $OUTPUT_DIR' ]; then
      out_file=${out_dir}/chunk_${chunk}.sh
      start_new_chunk ${out_file} ${header} ${STEP}

    elif [[ $line =~ ^STEP=.*$ ]]; then
      STEP=${line#STEP=}

    elif [[ $line =~ \#.JOB:.* ]]; then
      nb_job=$((nb_job+1))
      if [[ ${nb_job} -ge ${max_chunk} ]]; then
          nb_job=0  # reset counter
          echo '# END' >> $out_file
          chunk=$(($chunk+1))
          out_file=${out_dir}/chunk_${chunk}.sh
          start_new_chunk ${out_file} ${header} ${STEP}
          load_previous_submit_id ${out_file} ${out_dir}

          echo "$line" >> ${out_file}
      fi

    elif [[ $line =~ echo.*\ \>\>..JOB_LIST ]]; then
      job_list_chunk=${out_dir}/chunk_${chunk}.out
      bidon=$(echo "$line" | sed 's/echo\s"$\(.*\)\s$JOB_NAME.*/echo "\1=$\1 "/g')
      echo "$bidon >> ${job_list_chunk}" >> ${out_file}
    fi
done < ${genpipes_in}



# Submit and
