#!/bin/bash
#yanhong hong
#2020 - 12 - 09

energy_range=10


while getopts ":p:r:l:c:n:e:w:v:u:" OPT
do
  case $OPT in
    p) 
    process_num="$OPTARG";;
    r) 
    receptors_dir="$OPTARG";;
    l) 
    ligands_dir="$OPTARG";;
    c) 
    cpu_num="$OPTARG";;
    n) 
    num_modes="$OPTARG";;
    e) 
    exhaustiveness="$OPTARG";;
    w) 
    dockingdir="$OPTARG";;
    v) 
    vina_splitpath="$OPTARG";;
    u)
    uncomplete_tasks="$OPTARG";;    # 在本示例中，传进来的字符串大致如：
    ?) 
    echo "Invalid option -$OPTARG" >&2;;
  esac
done


ligandbasename=`basename ${ligands_dir}`
ligandpath=${ligands_dir}"/"${ligandbasename}".pdbqt"
sminapath="/home/qnhu/miniconda2/bin/smina"


mkdir ${dockingdir}/result
mkdir ${dockingdir}/log


for receptor in ${receptors_dir}/*.pdbqt
do
    receptorbasename=`basename $receptor .pdbqt`
    if [[ ${uncomplete_tasks} == *${receptorbasename}* ]]; then              # 在未完成任务字符串中提取出我们需要进行对接的任务受体id。
        rm -rf ${dockingdir}/result/${receptorbasename}_${ligandbasename}    # 强制删除非完整的任务。
        mkdir ${dockingdir}/result/${receptorbasename}_${ligandbasename}
        echo -e ${sminapath}" -r "${receptor}" -l "${ligandpath}" --autobox_ligand "${receptor}" --cpu "${cpu_num}" --num_modes "${num_modes}" --exhaustiveness "${exhaustiveness}" --energy_range "${energy_range}"  -o  "${dockingdir}/result/${receptorbasename}_${ligandbasename}/${receptorbasename}_${ligandbasename}.pdbqt" --log  "${dockingdir}/log/${receptorbasename}_${ligandbasename}.log  
    fi
done | xargs -n 1 -I {} -P ${process_num} bash -c '{}'

for result in ${dockingdir}/result/*/*.pdbqt
do
    resultbasename=`basename $result .pdbqt`
    receptorbasename=${resultbasename%_*}
    if [[ ${uncomplete_tasks} == *${receptorbasename}* ]]; then
      ${vina_splitpath} --input ${result}
    fi
done