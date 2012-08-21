###### Current path of the package 
cd ..
CPATH=$PWD
cd -

###### 0- Setting library ############
rm ~/.R/Makevars
#echo CXXFLAGS=-fpic -O3 -pipe -g -DOPTIMIZED > ~/.R/Makevars
R CMD INSTALL -l /home/vmiele/R/x86_64-pc-linux-gnu-library/2.14/ $CPATH/cghseg_1.0.0.tar.gz

####### 1- Parameters #################
M=50
N=100

####### 2- Creating test data #########
/usr/bin/time R --no-restore --no-save --args FALSE 1 simul $M $N TRUE FALSE <  cghseg.examples.args.R > /dev/null

####### 3- Integration test ###########
# with OPTIMIZATION (if enabled)
/usr/bin/time R --no-restore --no-save --args FALSE 1 load $M $N TRUE TRUE <  cghseg.examples.args.R | grep Integration | grep -v cat
# with PARRALEL 
/usr/bin/time R --no-restore --no-save --args FALSE 4 load $M $N TRUE FALSE <  cghseg.examples.args.R | grep Integration | grep -v cat
# with  PARRALEL+OPTIMIZATION (if enabled)
/usr/bin/time R --no-restore --no-save --args FALSE 4 load $M $N TRUE TRUE <  cghseg.examples.args.R | grep Integration | grep -v cat

##### 4- Profiling = turn the first ##
##### flag FALSE into TRUE     #######