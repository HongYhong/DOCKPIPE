######################################################################
#
#   a script to use plip to batch analysis protein-ligand interaction
#
#######################################################################

#author : yanhong hong



while getopts ":d:r:a:p:u:" OPT
do
  case $OPT in
    d) 
    resultdir="$OPTARG";;
    r) 
    receptordir="$OPTARG";;
    a) 
    analysisdir="$OPTARG";;
    p) 
    process_num="$OPTARG";;
    u)
    uncomplete_tasks="$OPTARG";;              #analysis_uncomplete_tasks 返回来的案例是 'complex_pdb1BHE_031a4cf825429ccfd5d341f50c40a039_ligand_01, complex_pdb1BHE_031a4cf825429ccfd5d341f50c40a039_ligand_02'
    ?) 
    echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done



if [ ! -d "$analysisdir" ]; then
    mkdir ${analysisdir}
fi


#处理文件组织格式，具体转化如下：
# '''
# |-- docking
# |   `-- 3.1.11.1_031a4cf825429ccfd5d341f50c40a039                                              
# |       |-- log
# |       |   |-- pdb1FXX_031a4cf825429ccfd5d341f50c40a039.log
# |       |   `-- pdb2QXF_031a4cf825429ccfd5d341f50c40a039.log          ----->             
# |       `-- result
# |           |-- pdb1FXX_031a4cf825429ccfd5d341f50c40a039.pdbqt
# |           `-- pdb2QXF_031a4cf825429ccfd5d341f50c40a039.pdbqt    
# '''

echo '##############################################'
echo 'analysis uncomplete tasks: '${uncomplete_tasks}
echo '##############################################'
    
# 把配体和蛋白结合起来
for result in ${resultdir}/result/*     
do
    result=`basename $result`
    receptorname=${result%_*} #extract the receptor name.In my case,e.g. A0A1D8AVC9_addH_DON --> A0A1D8AVC9
    if [[ $uncomplete_tasks == *"${result}"*  ]]; then
        echo '##########################################'
        echo 'analysis to do task:'${result}
        echo '##########################################'
        mkdir ${analysisdir}/${result}
        mkdir ${analysisdir}/${result}/receptor
        mkdir ${analysisdir}/${result}/complex
        cp ${receptordir}/${receptorname}.pdbqt ${analysisdir}/${result}/receptor
        cp ${resultdir}/result/${result}/*ligand_01.pdbqt ${analysisdir}/${result} # make sure that the result file name has "ligand"(e.g. A0A1D8AVC9_addH_DON_ligand_01.pdbqt) which is the default output of vina_split. 
        cd ${analysisdir}/${result}
        babel -ipdbqt receptor/${receptorname}*.pdbqt -opdb receptor/${receptorname}.pdb                                #总的流程：复制受体   --> 复制配体  --> 转化配体  --> 生成复合物 --> 创建和复合物同名的文件夹并进行分析。
        for ligand in *.pdbqt
        do
            ligandbasename=`basename $ligand .pdbqt`
            rm -rf ${ligandbasename}.pdb                #之前执行任务留下的pdb文件
            echo "babel -ipdbqt ${ligand} -opdb "${ligandbasename}".pdb"
        done | xargs -n 1 -I {} -P ${process_num} bash -c "{}"
    #combine receptor and ligand.
        for ligandpdb in *.pdb
        do
            ligandpdbbasename=`basename $ligandpdb .pdb`
            if [[ ${uncomplete_tasks} == *${ligandpdbbasename}* ]]; then
                rm -rf complex/complex_${ligandpdb}
                echo "babel -j "${ligandpdb}" receptor/"${receptorname}".pdb -opdb complex/complex_"${ligandpdb}   #耗时步骤
            fi
        done | xargs -n 1 -I {} -P ${process_num} bash -c "{}"
    fi
    
    
    #convert ligand file format.
done

#开始分析
for result in ${resultdir}/result/*
do
    result=`basename $result`
    if [[ ${uncomplete_tasks} == *${result}* ]]; then
        for complex in ${analysisdir}/${result}/complex/*.pdb
        do
            complex_base=`basename ${complex} .pdb`
            if [[ $uncomplete_tasks == *${complex_base}* ]]; then    
                rm -rf ${analysisdir}/${result}/complex/${complex_base}
                mkdir ${analysisdir}/${result}/complex/${complex_base}
                mv ${complex} ${analysisdir}/${result}/complex/${complex_base}
            fi
        done
    for complexdir in ${analysisdir}/${result}/complex/*
        do
            complexdirbasename=`basename ${complexdir}`
            if [[ $uncomplete_tasks == *${complexdirbasename}* ]]; then
                cd $complexdir
                echo -e ${complexdir}"\n docker run   -v "${PWD}":/results     -w /results     -u $(id -u "${USER}"):$(id -g "${USER}")  pharmai/plip:latest  -f *.pdb -xtyp"
                cd ${analysisdir}/..
            fi
        done | xargs -n 2  -P ${process_num}  -d'\n'  bash -c 'cd $0; $1'
    fi
done

#这个过程会产生哪些中间文件呢？  （1）${analysisdir}/${result} 下的 ligand 文件，只有一个 （2）受体文件 （3） ${analysisdir}/${result}/complex/${complex_base} 文件夹下的pdb文件
#analysis integrality 检验的是 txt 和 xml。

for result in ${resultdir}/result/*
do
    result=`basename $result`
    rm -rf ${analysisdir}/${result}/*.pdb*
    rm -rf ${analysisdir}/${result}/receptor/*.pdb*
done

