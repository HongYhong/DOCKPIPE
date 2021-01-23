# coding: utf-8

# author : yanhong hong

"""
只调回最佳构象的图片和Reprot
"""
import subprocess
from subprocess import PIPE, Popen
import os
from os.path import isdir, join
import time
import glob
from Crypto.Cipher import DES
from binascii import b2a_hex, a2b_hex
import hashlib
import zlib
import base64
import re
import shutil
import requests


CONVERT_ERROR = 1
DOCKING_ERROR = 2
ANALYSIS_ERROR = 3
GET_RESULT_ERROR = 4
COPYERROR = 5
PARAMERROR = 6




"""
FINISH 需要在前端展示的数据，如图片和pdb文件，存到 /home/qnhu/bsp/static/d_dock/results目录下

"""

"""
FINISH 删除 analysis 中间文件    在 analysis.sh 中执行
"""

def enzyme_recommand(ec, smiles, exhaustiveness=40, num_modes=10, docking_type='blinddock'):
    """
    如果已运行过该程序，则直接访问结果
    exhaustiveness，num_modes两个参数由用户自定义，将其加入文件夹命名规则中

    :param ec:
    :param smiles:
    :param exhaustiveness: 跟最后对接的精度有关系，越大准确度越高
    :param num_modes: 跟最后产生的构象数有关系
    :param docking_type:
    :return: 获取对接评分和图片
    """
    # ec = request.GET.get('ec')  # 酶类别
    # smiles = request.GET.get('smiles')  # 底物
    # xx = request.GET.get('x')  # 一些参数

    #  准备配体
    #  我们使用DES加密算法来加密我们的配体名称：
    key = b'SmilekeY'
    des = DES.new(key, DES.MODE_ECB)  # 创建一个DES实例
    smiles_hash = hashlib.md5(smiles).hexdigest()
    
    if ((num_modes > 99) or (num_modes < 10) or (num_modes % 1 != 0)):
        return PARAMERROR
    #plain_data = pad(str(smiles))  # 要加密的明文数据，长度必须是16的倍数
    
    #smiles_hash = base64.b32encode(smiles)
    #smiles_hash = smiles_hash.replace('=','~')

    workingdir = '/home/qnhu/bsp'
    # print(workingdir)
    workingdir = workingdir +'/d_dock'
    ligandsdir = workingdir + "/ligands"
    smilesdir = ligandsdir + "/" + str(smiles_hash)
    molsmi = smilesdir+ "/" + str(smiles_hash) + ".smi"
    molmol2 = smilesdir + "/" + str(smiles_hash) + ".mol2"
    molpdbqt = smilesdir+ "/" + str(smiles_hash) + ".pdbqt"

    # make directory for the input smile as a work directory for preparing ligands
    # 文件夹存在跳过转化过程
    smilesdir = smilesdir.strip()
    if not os.path.exists(smilesdir):
        print(smilesdir)
        os.mkdir(smilesdir)
        with open(molsmi, 'w') as molfile:
            molfile.write(smiles)

        # use babel to convert smile to mol2
        convertcmd = 'babel --gen3D -ismi {} -omol2 {}'.format(molsmi,molmol2)
        try:
            convert_status = subprocess.check_call(convertcmd, shell=True)
        except subprocess.CalledProcessError:
            return CONVERT_ERROR  # error code1

        #  mol2 to pdbqt
        pythonshpath = '/home/qnhu/MGLTools-1.5.6/mgltools_x86_64Linux2_1.5.6/bin/pythonsh'
        scriptpath = '/home/qnhu/MGLTools-1.5.6/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
        mol2topdbqt = pythonshpath + ' ' +scriptpath + ' -l ' + molmol2 + ' -o ' + molpdbqt

        try:
            convert_status2 = subprocess.check_call(mol2topdbqt, shell=True)
        except:
            return CONVERT_ERROR  # error code1

    # 开始对接

    #  配置参数：
    process_number = 1
    cpu_number = 1       # 总的cpu使用量 = process_number * cpu_number
    receptors_dir = workingdir + "/receptors/" + ec
    ligands_dir = workingdir + "/ligands/" + str(smiles_hash)
    docking_dir = workingdir + "/docking/" + ec + "_" + str(smiles_hash) + "_" + str(exhaustiveness) + "_" + str(num_modes)
    analysis_dir = workingdir + "/analysis/" + ec + "_" + str(smiles_hash) + "_" + str(exhaustiveness) + "_" + str(num_modes)
    vina_splitpath = '/home/qnhu/miniconda2/bin/vina_split'

    # 如果docking文件夹中已经存在我们运行过的任务那么直接返回结果。
    docking_uncomplete_tasks = dockingdir_integrality(docking_dir,receptors_dir,num_modes,ec,smiles_hash,exhaustiveness)
    print '########### docking uncomplete tasks ###############'
    print docking_uncomplete_tasks 
    print '##################################################'
    if len(docking_uncomplete_tasks) is not 0:    # FINISH 判断任务执行完整性函数
        if not os.path.exists(docking_dir):
            os.mkdir(docking_dir)

        #  强制删除 analysis dir 里面对应 uncomplete_tasks 中的目录。以免发生冲突。
        for uncomplete_task in docking_uncomplete_tasks:
            try:
                analysis_task_path = glob.glob(analysis_dir + '/*' + uncomplete_task + "*")[0]
                shutil.rmtree(analysis_task_path)
            except:
                pass

        print '######################### Analysis Result Files Removed #########################'


        #  对接命令
        #  判断对接模式：
        #  目前仅有盲对接
        docking_uncomplete_tasks = ','.join(map(str, docking_uncomplete_tasks))         #转换成字符串传入shell脚本。
        if docking_type == 'blinddock':  # 盲对接,即在整个受体中找到对接的位点。
            dockingcmd = 'bash /home/qnhu/bsp/d_dock/docking.sh -r "{}" -p{} -l "{}" -c{} -n{} -e{} -w "{}" -v "{}" -u "{}" '.format(
                receptors_dir, process_number, ligands_dir, cpu_number, num_modes, exhaustiveness, docking_dir, vina_splitpath, docking_uncomplete_tasks)
            f = open('/home/qnhu/bsp/d_dock/error_log.txt', 'w')
            f.write(str(dockingcmd))
            try:
                docking_status = subprocess.check_call(dockingcmd, shell=True)
            except subprocess.CalledProcessError:
                return DOCKING_ERROR  # error code2
        elif docking_type == 'topocket':  # 先找到靶点，再对接。
            pass
        else:
            pass

    analysis_uncomplete_tasks = analysisdir_integrality(docking_dir,receptors_dir,analysis_dir,smiles_hash)
    if len(analysis_uncomplete_tasks) is not 0:
        # analysis_uncomplete_task 返回来的案例是 complex_pdb1BHE_031a4cf825429ccfd5d341f50c40a039_ligand_01
        # FINISH 同上，判断任务执行完整性的函数。
        try:
            analysis_uncomplete_tasks = ','.join(map(str, analysis_uncomplete_tasks))
            analysis_result(receptors_dir, docking_dir, analysis_dir, process_number, analysis_uncomplete_tasks)
        except:
            # docking 文件夹中存在， analysis 文件夹中不存在
            return ANALYSIS_ERROR

    
    # 得到 pdb评分，Report，和相互作用图
    print(analysis_dir)

    # try:
    result_list = get_results(analysis_dir,docking_dir)
    # except:
    #     return GET_RESULT_ERROR
    return result_list


def pad(text):
    """加密函数, 加密文本text必须为8的倍数
    如果text不是8的倍数，那就补足为8的倍数"""
    while len(text) % 8 != 0:
        text += ' '
    return text


def analysis_result(receptors_dir, docking_dir, analysis_dir, process_number, analysis_uncomplete_tasks):
    analysiscmd = 'bash /home/qnhu/bsp/d_dock/analysis.sh -d "{}" -r "{}" -a "{}" -p{} -u "{}"'.format(
        docking_dir, receptors_dir, analysis_dir, process_number, analysis_uncomplete_tasks)
    try:
        analysis_status = subprocess.check_call(analysiscmd, shell=True)
    except subprocess.CalledProcessError:
        return ANALYSIS_ERROR


def get_score(parentdir,poserank,docking_dir):
    pdbid = parentdir.split('_')[1]
    smiles_hash = parentdir.split('_')[2]
    poserank = re.sub("^0","",poserank)
    logpath = docking_dir + "/log/" + pdbid + "_" + smiles_hash + ".log"
    get_scorecmd = 'grep "' + '^' + poserank + ' " ' + logpath + "|awk '{print $2}'"
    try:
        score = subprocess.check_output(get_scorecmd,shell=True).strip('\n')
    except subprocess.CalledProcessError:
        return GET_SCORE_ERROR
    return score


def get_results(analysis_dir,docking_dir):
    analysispaths = [join(analysis_dir, d) for d in os.listdir(analysis_dir) if isdir(join(analysis_dir, d))]
    print(analysispaths)
    # analysis 文件夹中需要返回的包括 png, pse, report.txt, report.xml
    reporttextpaths = []
    reportxmlpaths = []
    ec = os.path.basename(analysis_dir.split('_')[1])


    print 'analysisdir',analysis_dir
    print 'in get results############################################################ ',ec


    ##############################################################################重命名部分#################################################################################

    # 返回所有 report.txt 的绝对路径
    for dirpath, subdirs, files in os.walk(analysis_dir):
        for x in files:
            if x.endswith(".txt"):
                reporttextpaths.append(os.path.join(dirpath, x))    
                    
    # 由于标准输出的 report.txt的命名不具有特异性，我们需要对这些 txt进行重命名，使用
    # report.txt 的上一级目录名称。
    for reporttextpath in reporttextpaths:
        reportdir = os.path.abspath(os.path.join(reporttextpath, os.pardir))
        # e.g. /home/qnhu/bsp/yanhong/analysis/3.2.1.15_031a4cf825429ccfd5d341f50c40a039/pdb1K5C_031a4cf825429ccfd5d341f50c40a039/complex/complex_pdb1K5C_031a4cf825429ccfd5d341f50c40a039_ligand_02

        parentdir = os.path.basename(os.path.normpath(reportdir))
        os.rename(reporttextpath,reportdir + "/" + parentdir + "_report.txt")

    reporttextpaths = []

    # 重新赋值
    for dirpath, subdirs, files in os.walk(analysis_dir):
        for x in files:
            if x.endswith(".txt"):
                reporttextpaths.append(os.path.join(dirpath, x))

    # ###########################################################################
    # 获取所有xml文件
    for dirpath, subdirs, files in os.walk(analysis_dir):
        for x in files:
            if x.endswith(".xml"):
                reportxmlpaths.append(os.path.join(dirpath, x)) 
    
    for reportxmlpath in reportxmlpaths:
        reportdir = os.path.abspath(os.path.join(reportxmlpath, os.pardir)) #e.g. /home/qnhu/bsp/yanhong/analysis/3.2.1.15_031a4cf825429ccfd5d341f50c40a039/pdb1K5C_031a4cf825429ccfd5d341f50c40a039/complex/complex_pdb1K5C_031a4cf825429ccfd5d341f50c40a039_ligand_02
        parentdir = os.path.basename(os.path.normpath(reportdir))
        os.rename(reportxmlpath, reportdir + "/" + parentdir + "_report.xml")

    reportxmlpaths = []

    for dirpath, subdirs, files in os.walk(analysis_dir):
        for x in files:
            if x.endswith(".xml"):
                reportxmlpaths.append(os.path.join(dirpath, x))


    ########################################重命名部分结束############################################################################

    result_list = []

    for reportxmlpath in reportxmlpaths:
        try:

            origin_dir =  '/home/qnhu/bsp/d_dock/analysis/'
            static_dir = '/home/qnhu/bsp/static/d_dock/results/'
            reportdir = os.path.abspath(os.path.join(reportxmlpath, os.pardir))
            parentdir = os.path.basename(os.path.normpath(reportdir))

            pdbid = parentdir.split('_')[1]
            ecid = ec
            pose_rank = parentdir.split('_')[-1]
            imgpath = glob.glob(reportdir + "/*_LIG_*.png")[0]  #e.g. imgpath: /home/qnhu/bsp/yanhong/analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/pdb2QXF_031a4cf825429ccfd5d341f50c40a039/complex/complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_05/COMPLEX_PDB2QXF_031A4CF825429CCFD5D341F50C40A039_LIGAND_05_PROTEIN_LIG_d_1.png
            imgstaticpath = imgpath.replace(origin_dir,static_dir) #我们把它拷贝到：imgstaticpath:  /home/qnhu/bsp/static/d_dock/results/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/pdb2QXF_031a4cf825429ccfd5d341f50c40a039/complex/complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_05/COMPLEX_PDB2QXF_031A4CF825429CCFD5D341F50C40A039_LIGAND_05_PROTEIN_LIG_d_1.png
            copyimgcmd = 'mkdir -p ' + os.path.dirname(imgstaticpath) + ' && cp ' + imgpath + ' ' + os.path.dirname(imgstaticpath)
            try:
                copyimgstatus = subprocess.check_call(copyimgcmd, shell=True)
            except:
                print(COPYERROR)
                # return COPYERROR
            psepath = glob.glob(reportdir + "/*_LIG_*.pse")[0] # pse : pymol session file.
            psestaticpath = psepath.replace(origin_dir,static_dir)
            copypsecmd = 'mkdir -p ' + os.path.dirname(psestaticpath) + ' && cp ' + psepath + ' ' + os.path.dirname(psestaticpath)
            try:
                copypsestatus = subprocess.check_call(copypsecmd, shell=True)
            except:
                print(COPYERROR)
                # return COPYERROR

            pdbpath = glob.glob(reportdir + "/*protonated.pdb")[0]
            pdbstaticpath = pdbpath.replace(origin_dir,static_dir)
            copypdbcmd = 'mkdir -p ' + os.path.dirname(pdbstaticpath) + ' && cp ' + pdbpath + ' ' + os.path.dirname(pdbstaticpath)
            try:
                copypdbstatus = subprocess.check_call(copypdbcmd, shell=True)
            except:
                print(COPYERROR)


            reporttextpath = glob.glob(reportdir + "/*.txt")[0]
            reporttextstaticpath = reporttextpath.replace(origin_dir,static_dir)
            copyreporttextcmd = 'mkdir -p ' + os.path.dirname(reporttextpath) + ' && cp ' + reporttextpath + ' ' + os.path.dirname(reporttextstaticpath)
            try:
                copyreporttextstatus = subprocess.check_call(copyreporttextcmd, shell=True)
            except:
                print(COPYERROR)
                # return COPYERROR
            reportxmlpath = reportxmlpath
            reportxmlstaticpath = reportxmlpath.replace(origin_dir,static_dir)
            copyreportxmlcmd = 'mkdir -p ' + os.path.dirname(reportxmlpath) + ' && cp ' + reportxmlpath + ' ' + os.path.dirname(reportxmlstaticpath)
            try:
                copyreportxmlstatus = subprocess.check_call(copyreportxmlcmd,shell=True)
            except:
                print(COPYERROR)
                # return COPYERROR
            

            print(parentdir)
            print(pose_rank)
            print(docking_dir)
            print "i am here #####################################################################################1"
            score = get_score(parentdir, pose_rank, docking_dir)
            print "i am here #####################################################################################2"
            print(score)


            print "i am here #####################################################################################3"
            print pdbstaticpath
            result_dct = {
                "pdbid" : pdbid,
                "psepath": psepath,
                "ec": ecid,
                "pose_rank": pose_rank,
                "pdbpath": pdbstaticpath,
                "imgpath": imgstaticpath,
                "psepath": psestaticpath,
                "reporttextpath": reporttextstaticpath,
                "reportxmlpath": reportxmlstaticpath,
                "score": score,
            }

            result_list.append(result_dct)
        except:
            print("Error")

    return result_list

# 检查已完成任务完整性
def dockingdir_integrality(docking_dir,receptors_dir,num_modes,ec,smiles_hash,exhaustiveness):
    receptors_ids = [ receptor.replace(".pdbqt","") for receptor in os.listdir(receptors_dir) if os.path.isfile(join(receptors_dir ,receptor))]
    logspath_list = [ docking_dir + "/log/" + receptor + "_" + smiles_hash + ".log" for receptor in receptors_ids ]
    resultspath_list = [ docking_dir + "/result/" + receptor + "_" + smiles_hash + ".pdbqt" for receptor in receptors_ids ]
    
    docking_uncomplete_tasks = set()
    #check log files existance.检查log文件是否存在
    for log in logspath_list:
        found = False
        if not os.path.exists(log):
            logbasename = os.path.basename(log).split('_')[0]
            docking_uncomplete_tasks.add(logbasename)   # FINISH: 非完整的文件结构仅返回执行未完毕的任务set。 UNCOMPLETE_TASK = {}
        else:       # chech the integrality of the log files.
            with open(log,'r') as logfile:
                for line in logfile:
                    if re.search("^" + str(num_modes) ,line):
                        found = True
                        break
            if found == False:
                logbasename = os.path.basename(log).split('_')[0]
                docking_uncomplete_tasks.add(logbasename)

    #check result files existance.
    for result in resultspath_list:
        found = False
        if not os.path.exists(result):
            resultbasename = os.path.basename(log).split('_')[0]
            docking_uncomplete_tasks.add(resultbasename)
        else:       # check integrality of the result files.
            with open(result,'r') as resultfile:
                for line in resultfile:
                    if re.search("^MODEL " + str(num_modes),line):
                        found = True
                        break
            if found == False:
                resultbasename = os.path.basename(log).split('_')[0]
                docking_uncomplete_tasks.add(resultbasename)

    return docking_uncomplete_tasks

# 20210115 设定 num_modes = 1 仅展示一个结果即可.

def analysisdir_integrality(docking_dir,receptors_dir,analysis_dir,smiles_hash,num_modes = 1):  # FINISH 函数测试完毕 ，测试模块 # test for analysisdir integrality.
    analysis_uncomplete_tasks = []
    reporttextfilelist = []
    reportxmlfilelist = []

    receptors_ids = [ receptor.replace(".pdbqt","") for receptor in os.listdir(receptors_dir) if os.path.isfile(join(receptors_dir,receptor)) ]
    print '############# print receptorids #####################',receptors_ids
    # print("#############print complex_dir_prefixs##############\n",complex_dir_prefixs)
    print "params:",analysis_dir,receptor,smiles_hash
    complex_dir_prefixs = [ analysis_dir + "/" + receptor + "_" + smiles_hash + "/complex/" + "complex_" + receptor + "_" + smiles_hash + "_ligand_" for receptor in receptors_ids ]
    
    
    complex_dirs = []
    for complex_dir_prefix in complex_dir_prefixs:
        for index in range(num_modes):
            suffix = "%02d" % (index+1,)
            complex_dir = complex_dir_prefix + suffix
            complex_dirs.append(complex_dir)
    #check the existances of pymol session file, binding site view png file,interaction report text file, interaction report xml file. # ! 注意 plip 的分析如果没有相互作用，那么
                                                                                                                                       # ! 输出结果是没有 pse 文件和 png 文件的，所以判断分析是否完整我们只看txt文件和xml文件。
    for complex_dir in complex_dirs:  #这个检测的前提是如果检测到 txt 文件和 xml 文件，不管其内容是否完整，我们都认为它是一个完成的任务。
        try:
            reporttextfilelist.append(glob.glob(complex_dir + "/*.txt")[0])
            reportxmlfilelist.append(glob.glob(complex_dir + "/*.xml")[0])
        except:
            analysis_uncomplete_tasks.append(os.path.basename(complex_dir))

    return analysis_uncomplete_tasks



if __name__ == "__main__":
    # tic = time.time()

    # ec = '3.2.1.15'
    # # ec = '3.1.11.1'
    # smiles = 'C1=CC=C(C=C1)C=O'
    # results = enzyme_recommand(ec, smiles)
    # #
    # toc = time.time()
    # print(results)
    # print(toc - tic)

    # test for analysisdir integrality.
    # analysis_uncomplete_task  = analysisdir_integrality('/home/qnhu/bsp/yanhong/docking/3.2.1.15_031a4cf825429ccfd5d341f50c40a039','/home/qnhu/bsp/yanhong/receptors/3.2.1.15',40,'/home/qnhu/bsp/yanhong/analysis/3.2.1.15_031a4cf825429ccfd5d341f50c40a039','031a4cf825429ccfd5d341f50c40a039')
    # print(analysis_uncomplete_task)

    # test for dockingdir integrality.
    # param docking_dir,receptors_dir,num_modes,analysis_dir,smiles_hash
    #docking_uncomplete_tasks =  dockingdir_integrality('/home/qnhu/bsp/yanhong/docking/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40','/home/qnhu/bsp/yanhong/receptors/3.1.11.1',40,'/home/qnhu/bsp/yanhong/analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039','031a4cf825429ccfd5d341f50c40a039',40)
    # docking_uncomplete_tasks =  dockingdir_integrality('/home/qnhu/bsp/yanhong/docking/3.2.1.15_031a4cf825429ccfd5d341f50c40a039','/home/qnhu/bsp/yanhong/receptors/3.2.1.15',40,'/home/qnhu/bsp/yanhong/analysis/3.2.1.15_031a4cf825429ccfd5d341f50c40a039','031a4cf825429ccfd5d341f50c40a039',40)
    # print docking_uncomplete_tasks

    #  开始调试对接任务的执行。 2020/12/18
    # FINISH 对接脚本 docking.sh 分析脚本 analysis.sh 的修改。

    # FINISH get_results 函数的调试。
    # test for get_results:
    # param analysis_dir, docking_dir
    # results_list = get_results("/home/qnhu/bsp/yanhong/analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40","/home/qnhu/bsp/yanhong/docking/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40")
    # print(results_list)

    # FINISH analysis_result 函数调试
    # test for analysis_results:
    # param receptors_dir, docking_dir, analysis_dir, process_number, analysis_uncomplete_tasks
    #analysis_result("/home/qnhu/bsp/yanhong/receptors/3.1.11.1", "/home/qnhu/bsp/yanhong/docking/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40", "/home/qnhu/bsp/yanhong/analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40", 1, "complex_pdb1FXX_031a4cf825429ccfd5d341f50c40a039_ligand_01, complex_pdb1FXX_031a4cf825429ccfd5d341f50c40a039_ligand_02, complex_pdb1FXX_031a4cf825429ccfd5d341f50c40a039_ligand_03 , complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_01 ,complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_02,complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_03,complex_pdb2QXF_031a4cf825429ccfd5d341f50c40a039_ligand_05")
    

    # TODO 总体调试
    # FINISH 测试成功 test1 : 测试非完整任务的执行情况：  
    '''
        before:
        tree analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/ > testtext/beforetest1analysis.txt
        tree docking/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/ > testtext/beforetest1docking.txt
    '''
    # enzyme_recommand('3.1.11.1', 'C1=CC=C(C=C1)C=O', exhaustiveness=40, num_modes=40, docking_type='blinddock')

    '''
        after:
        tree analysis/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/ > testtext/aftertest1analysis.txt
        tree docking/3.1.11.1_031a4cf825429ccfd5d341f50c40a039_40_40/ > testtext/aftertest1docking.txt
    '''
    
    # FINISH test2 : 测试完全新任务的执行情况：
    #ec = '3.1.11.1'
    #smiles = 'C1=CC=C(C=C1)C=O'
    #results_list = enzyme_recommand(ec, smiles, exhaustiveness=30, num_modes=10, docking_type='blinddock')
    #print(results_list)
    #失败  文件名过长！ OSError: [Errno 36] File name too long: './plipfixed.complex_pdb3C94_3c37262d63368b73bcf4a391dd66782e6a10c5d9a5f453dec28541893aa6963006413dd094a4b3e8216b578b54a7d35de83975d29399a5448aef6ca8f68bad826334eaea878cc5b3ab92d027234e9f1e67c4074be1cdd6c6a86103a529b9624ca2b87693e552c923_ligand_01_xwsw76fw.pdb'
    #修改： 使用其它的加密方式 使得输出长度是一定的。   尝试使用 zlib 压缩 + base64 进行 encoding. # ! 此时不具有加密功能,仅仅是压缩字符串。 
    # 已修复

    #FINISH test3: 测试完全新任务执行情况2：
    # ec = '3.2.1.15'
    # smiles = 'CC1(C)O[C@H]2C[C@@H]3[C@@H]4C[C@H](F)C5=CC(=O)C=C[C@]5(C)[C@]4(F)[C@H](O)C[C@@]3(C)[C@]2(C(=O)CO)O1'
    # results_list = enzyme_recommand(ec, smiles, exhaustiveness = 30 ,num_modes = 10 ,docking_type= 'blinddock')
    # print(results_list)

    # 有 bug，因为如果原受体中含有其它的配体，那么plip输出结果会有多个png与pse，我们要仅返回我们输入的配体和受体之间的相互作用的png 和 pse。
    # bug已修复

    # FINISH test4 : 测试完成任务的执行情况：
    # ec = '3.2.1.15'
    # smiles = 'CC1(C)O[C@H]2C[C@@H]3[C@@H]4C[C@H](F)C5=CC(=O)C=C[C@]5(C)[C@]4(F)[C@H](O)C[C@@]3(C)[C@]2(C(=O)CO)O1'
    # results_list = enzyme_recommand(ec, smiles, exhaustiveness = 30 ,num_modes = 10 ,docking_type= 'blinddock')
    # print(results_list)

    # TODO test5:

    ec = '3.2.1.15'
    smiles = 'C[C@@]1(c2ccccc2)OC(C(=O)O)=CC1=O'
    results_list = enzyme_recommand(ec ,smiles, exhaustiveness = 10 ,num_modes = 10 ,docking_type = 'blinddock')
    print(results_list)