---
title: How to develop a Python script for DFT calculation job
description: This is from my jianshu platform. 
category: Sample
image: A_ship.jpg
sitemap: True
---

在上adv. metal课的时候就曾经遇到过这个问题，当然那个时候还是全靠手动去做nuclei move的锅，不过现在也会遇到这样的情况，就是我要提交一系列的任务，而这一系列的任务是有一定的规律的，目前从我遇到过的case来看，主要有两种

# Case 1. 
这种情况每个任务里面不一样的地方只有其中一个参数，参数是单一的或者多个参数之间有某种数学关系
e.g. 要做一个convergence test，提交相同cut off energy不同k-point的一系列任务，或者POSCAR里面成比例地修改键长

用改POSCAR这个case来说，需要完成两个script：
1. 生成POSCAR的script
首先理清思路，在我们的一系列的任务中，POSCAR是一直在变的，所以在最开始的输入文件中，我们不事先提供POSCAR这个文件，而是用每个文件夹中的script去生成跟任务文件夹对应的POSCAR文件，然后运行vasp程序->得到结果->收集结果

2. 写一个生成POSCAR的script pos.py，这里用python写了一个，用c也可以，完全看习惯

```
import sys
import math

a = float(sys.argv[1]) # lattice constant
fout = open('POSCAR', 'w') # make a file named as POSCAR

fout.write(
'''G_monolayer
1.0
        %.10f         0.0000000000         0.0000000000
        %.10f         %.10f         0.0000000000
        0.0000000000         0.0000000000         3.0000000000
    C
    2
Cartesian
     0.000000000         0.000000000         2.171259403
     1.233863017         0.712369702         2.171259403
'''%(a, -a/2, math.sqrt(3)*a/2))

'''construct a primitive cell with constant c direction but variable a b directions'''

fout.close()
```

这里我们是为了改变一个graphene的primitive cell的lattice parameter，由于是一个固定的Triclinic的结构，所以a和b的关系是固定的，也就是说只相当于一个变量，每给出一个a值，我们就可以用这个script得到一个对应的POSCAR

3. 提交job的script
在完成了POSCAR的generator之后，写另外一个script job用来调用前面的script去生成一系列的POSCAR并运行vasp程序

```
#!/bin/bash

#PBS -N GRAPHENE
#PBS -q USERNAME 
#PBS -l nodes=2:ppn=24
#PBS -l walltime=48:00:00
#PBS -V
#PBS -S /bin/bash

cd $PBS_O_WORKDIR

#source /public/software/profile.d/openmpi-intel-env.sh
#source /public/software/profile.d/intel-env.sh

source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/intel64/mklvars_intel64.sh
source /opt/intel/impi/5.0.2.044/bin64/mpivars.sh

EXEC=/opt/software/vasp/vasp-5.3.5-base
#EXEC=/opt/software/vasp/vasp-5.3.5-noZ

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | wc -l`

#partition line, from here start to call previous script and run vasp program

rm WAVECAR SUMMARY
for i in $(seq 2.4590 0.0005 2.4690)
do
     python pos.py $i

     # relaxation
     mpirun -genv I_MPI_DEVICE rdma -machinefile $PBS_NODEFILE -n $NP $EXEC > vasp.log

     E=`awk '/F=/ {print $0}' OSZICAR` ; echo $i $E >>SUMMARY
     mkdir $i
     mv CHG CONTCAR EIGENVAL OSZICAR PCDAT vasprun.xml WAVECAR CHGCAR DOSCAR IBZKPT OUTCAR POSCAR XDATCAR vasp.log $i

done
```

这个script的工作流程就是，对于一个特定的i值->调用之前的pos.py生成POSCAR->用生成的POSCAR搭配之前提供的INCAR等文件进行vasp计算，并且将计算结果（一个单步的计算，所以只有一个离子步的结果）放入名为SUMMARY的文件中，然后把计算过程中输出的文件和当前的POSCAR一起放入命名为值的文件夹中

运行之前要提交的文件:

![输入的文件](https://upload-images.jianshu.io/upload_images/20672840-3daa412ef041c6e9.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/720)

qsub job 运行之后应该会得到这样的结果：

![输出的文件](https://upload-images.jianshu.io/upload_images/20672840-3d48f6f1eeb70175.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/720)

然而实际上我们只需要最后下载SUMMARY这个文件就够啦

```
2.4590 1 F= -.16200431E+02 E0= -.16200431E+02 d E =-.377889E-09
2.4595 1 F= -.16200510E+02 E0= -.16200510E+02 d E =-.160482E-08
2.4600 1 F= -.16200576E+02 E0= -.16200576E+02 d E =-.623061E-08
2.4605 1 F= -.16200619E+02 E0= -.16200619E+02 d E =-.221167E-07
2.4610 1 F= -.16200653E+02 E0= -.16200653E+02 d E =-.717880E-07
2.4615 1 F= -.16200667E+02 E0= -.16200667E+02 d E =-.213098E-06
2.4620 1 F= -.16200671E+02 E0= -.16200671E+02 d E =-.578529E-06
2.4625 1 F= -.16200661E+02 E0= -.16200661E+02 d E =-.143590E-05
2.4630 1 F= -.16200629E+02 E0= -.16200627E+02 d E =-.325742E-05
2.4635 1 F= -.16200580E+02 E0= -.16200576E+02 d E =-.674955E-05
2.4640 1 F= -.16200515E+02 E0= -.16200509E+02 d E =-.127662E-04
2.4645 1 F= -.16200445E+02 E0= -.16200434E+02 d E =-.220353E-04
2.4650 1 F= -.16200355E+02 E0= -.16200338E+02 d E =-.347476E-04
2.4655 1 F= -.16200252E+02 E0= -.16200227E+02 d E =-.501850E-04
2.4660 1 F= -.16200157E+02 E0= -.16200124E+02 d E =-.666691E-04
2.4665 1 F= -.16200046E+02 E0= -.16200005E+02 d E =-.818921E-04
2.4670 1 F= -.16199922E+02 E0= -.16199875E+02 d E =-.934967E-04
2.4675 1 F= -.16199788E+02 E0= -.16199739E+02 d E =-.996283E-04
2.4680 1 F= -.16199656E+02 E0= -.16199606E+02 d E =-.992796E-04
2.4685 1 F= -.16199505E+02 E0= -.16199459E+02 d E =-.925687E-04
2.4690 1 F= -.16199381E+02 E0= -.16199341E+02 d E =-.805424E-04
```

然后复制到excel里面进行处理就好了，当然这里之所以可以这样的原因主要还是因为这个计算只有一个离子步，所以每一次运算只有一个E0=...的输出，如果有多个离子步的计算就需要作相应的调整从而只抓取最后一个E0的值了


这种情况感觉只适用于那些变量比较单一，而且有一定的数学规律的情况，遇到那种结构优化过之后的复杂结构，输入的POSCAR都无规律的不一样的情况就比较无力了

# Case 2.
对于比较复杂的结构，可以先在本地计算机上建立好文件夹，然后用script单纯地批量提交任务而不去改变文件里面内容，当然用这种方式得到的结果也只是一个个的独立的文件夹，获取结果的时候需要另外的script或者直接手动-_-||
用优化过的结构的bader电荷分析为例：
1. 复制之前的含有结果的文件夹，mv CONTCAR POSCAR，然后修改INCAR的参数到适合bader电荷分析的计算，上传到服务器上，这里的文件是12-20的一组数字命名的文件夹

2. 登陆服务器，在包含一批任务的父文件夹中touch jobsub，然后vi jobsub编写这个script

```
#!/bin/bash
for i in $(seq 12 1 20)
do
cd $i/
qsub job-gamma
cd ../
done
```

然后进入包含一系列任务的父文件夹中

![包含一系列任务文件夹的父文件夹](https://upload-images.jianshu.io/upload_images/20672840-846e8670d0f58323.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/720)

运行这个script，这里简单说一下bash的应用，这样写出来的一个script应该是".sh"的一个文件，也就是说第一行写的东西为这个文件标明了身份（用什么语言去compile并执行），第一行写的是bash，于是就是用bash来运行，如果单独运行这个文件的话，就像这样

```
[xxx@xxx double]$ ls
12  13  14  15  16  17  18  19  20  jobsub
[xxx@xxx double]$ bash jobsub 
271306.mu01
271307.mu01
271308.mu01
271309.mu01
271310.mu01
271311.mu01
271312.mu01
271313.mu01
271314.mu01

```

就好了
这个script原理十分简单，里面的代码就是一个loop，其实就是用这个操作代替了我们手工的操作，节约了大量的时间，尤其是这种用手工作容易出错且浪费时间的工作