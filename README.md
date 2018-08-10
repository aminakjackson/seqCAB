# seqCAB
Sequence Convertor &amp; Annotation with BLAST+

seqCAB is a BLAST utility to optimize high-throughput ‘omics’ analyses through multi-core parallelization and taxonomic annotation" Jackson et al. (In review). It uses MPI to allocate resources for increased running efficiency on both single and multiple nodes. 

This repository hosts a source code for seqCAB MPI and seqCAB annotation python scripts, FASTA protein and nucleotide test datasets (described in Test Dataset description). It has been tested on SGI and CRAY HPC systems, Ubuntu and Mac (only annotation).  In depth usage and running of seqCAB MPI and annotation can be found in pdf file is also included in this repository under documentation.  
 
#Getting seqCAB
The most up-to-date versions of seqCAB MPI and seqCAB annotation can be downloaded from out Git repository. To accomplish this, you can execute on your shell. 
$ git clone https://github.com/aminakjackson/seqCAB.git  
For updates, enter seqCAB directory and execute:
$ git pull
You can also click download zip to have seqCAB in your download directory. 


#############################
## Compiling BLAST+ 
#############################

1. Download NCBI BLAST from NCBI site
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-src.tar.gz

2. Uncompress the tar file:
$ tar xvzf ncbi-blast-2.2.31+-src.tar.gz
Enter the created NCBI directory:
$ cd ncbi-blast-2.2.31+-src

3. Reveal the compilers in your environment 
$ echo $CXX
$ echo $CC

If these are not discoverable, you may need to install Intel® Parallel Studio XE. It could also be that the modules with these compliers are not loaded.  

Note: Please do not compile your BLAST with a gcc compiler. MPI does not run with these. 
We noticed that the latest BLAST+ versions do not seem to compile with icc or mpiicc compliers. Our tests worked with blast+/2.2.31/ncbi-blast-2.2.31+-src and blast+/2.2.31/ncbi-blast-2.2.30+-src

When compiling, you have two options. You can either use mpiicc & mpiicpc or icc & icpc 
Depending on your system configurations: 
For SGI system you may need to load the intel-mpi complier module if not loaded. If loaded skip to step 7, if unsure perform steps 4 -6.
4. To see what modules are loaded, perform
$ module list
If the intel-mpi module is loaded skip to step 7 if not loaded:
5. To see what intel-mpi modules are available execute
$ module avail 
6. To load the module execute
$ module load intel-mpi**

If running on Cray system, the default is usually PrgEnv-cray and c++ compliers CC and cc. BLAST could not compile with these. You need to perform step 4 and copy your PrgEnv-cray module. Execute 5 to see what intel modules you have available. Then instead of 6 execute
$ module swap PrgEnv-cray PrgEnv-intel
Check if the icc complier is accessible by executing 
$ which icc
if not accessible execute 5 again and load an intel-mpi module that is available.

7. If you are planning on running on Xeon Host and Xeon Phi, you will need to compile for both. It is recommended that you duplicate the c++ directory before compiling for simplicity. Execute 
$ cp -rp c++ icc or cp -rp mpiicc

8. Change the directory <ncbi-blast-2.2.31+-src/icc$ or <ncbi-blast-2.2.31+-src/mpiicc$ and run
$ ./configure --with-bin-release --without-debug --without-gui --with-mt --without-boost --without-z --without-bz2 --without-strip CXX=icpc CC=icc --with-symbols

Or 
$ ./configure --with-bin-release --without-debug --without-gui --with-mt --without-boost --without-z --without-bz2 --without-strip CXX=mpiicpc CC=mpiicc --with-symbols

If your configuration is unsuccessful BLAST usually has issues with lower gcc versions. Check you gcc complier
$ which gcc. 
If the version the lower than 5 perform execute, step 5 to see which ones are available and then swap it to a version higher than 5 if available by running.
$ module swap old_gcc new_gcc then repeat 8

9. Modify src/algo/blast/core/blast_gapalign.c by inserting a line before function  definition BLAST_GetGappedScore  at code line 3401 in BLAST version 31


10. Change directory to ReleaseMT/build
$ vi Makefile.mk and edit by Changing the -O2 flag on lines 93 and 94 to -O3, and on 98 change -O to -O3.

11. Still in the build directory run
$ make all_r or make -k all_r -j8 

make -k all_r -j8 is much faster but it ignores the compilation errors that you will have to check the config file in case of unsuccessful compilation 

For successful build, check bin directory to see if all the BLAST executables are there.


Compiling for Xeon Phi

1. You will need to request for an interactive session on the Phi nodes. Repeat steps 3 to 6.

2. Run 
$. /configure --with-bin-release --without-debug --without-gui --with-mt --without-boost --without-z --without-bz2 --without-strip CXX=icpc CC=icc --with-symbols --with-build-root=ReleasePHI

Or 
./configure --with-bin-release --without-debug --without-gui --with-mt --without-boost --without-z --without-bz2 --without-strip CXX=mpiicpc CC=mpiicc --with-symbols --with-build-root=ReleasePHI

3. Change directory to ReleaseMT/build
$vi Makefile.mk and edit by Changing the -O2 flag on lines 93 and 94 to -O3 -mmic, and on 98 change -O to -O3 -mmic.

Provide the path to NCBI data tools by exporting it in the environment 
$ export NCBI_DATATOOL_PATH=<PATH$/ncbi-blast-2.2.30+-src/icc/ReleaseMT/bin or export NCBI_DATATOOL_PATH=/ncbi-blast-2.2.30+-src/mpiicc/ReleaseMT/bin

Make the following modifications in /ncbi-blast-2.2.31+-src/c++/include/ and edit corelib/ncbifloat.h  on line 71
#  if __cplusplus $= 201103L && defined(_GLIBCXX_CONSTEXPR)
Add 
&&! defined (__MIC__)



4. Change directory to ReleasePHI/build and run
$ make all_r or make -k all_r -j8 

Look in bin for executables.

#################################
## Installing newer versions of anaconda
##################################

 Downloading (Anaconda3 is recommended although 2 would work as well). 
1.	$ wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh 

#bash Anaconda3-5.2.0-Linux-x86_64.sh ## takes about 10 mins for all packages to install
## need to start a new terminal after installation is complete
#
-	Installer prompts “In order to continue the installation process, please review the license agreement.” Click Enter
-	Scroll to the bottom of the license terms and enter “Yes” to agree.
-	The installer prompts you to click Enter to accept the default install location, CTRL-C to cancel the installation, or specify an alternate installation directory. If you accept the default install location, the installer displays “PREFIX=/home/<user>/anaconda<2 or 3>” and continues the installation. It may take a few minutes to complete.
-	 The installer prompts “Do you wish the installer to prepend the Anaconda<2 or 3> install location to PATH in your /home/<user>/.bashrc ?” Enter Yes.

Testing Anaconda installation:
There are various ways to test successful anaconda installation after completion 
1.	Opening Anaconda Navigator, a program that is included with Anaconda by typing anaconda-navigator in a new shell.  
Note: If your system does not support graph interface, this may not work. 
2.	You can also try cd /anaconda3/bin and executing conda update conda
Note: without sudo permissions, Anconda is likely installed in /home/user/anaconda3/bin.
If step 2 fails try $ source activate conda this may display an error message but it would still do the trick
3.	To see what packages are installed execute 
$ conda list
In newer versions of Anaconda, pandas package is preinstalled. 
2. Use conda to install other packages
- $ conda install mpi4py
-  $conda install biopython

