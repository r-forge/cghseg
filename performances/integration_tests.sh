###### Current path of the package 
CPATH=$PWD

###### Setting library ############
rm ~/.R/Makevars
#echo CXXFLAGS=-fpic -O3 -pipe -g -DOPTIMIZED > ~/.R/Makevars
R CMD INSTALL -l /home/vmiele/R/x86_64-pc-linux-gnu-library/2.14/ $CPATH/cghseg_0.0.1.tar.gz

####### Parameters #################
M=50
N=100
cd $CPATH/performances

####### Creating test data #########
/usr/bin/time R --no-restore --no-save --args FALSE FALSE simul $M $N TRUE FALSE <  cghseg.examples.args.R > /dev/null

###### Integration test ############
# with OPTIMIZATION
/usr/bin/time R --no-restore --no-save --args FALSE FALSE load $M $N TRUE TRUE <  cghseg.examples.args.R | grep Integration | grep -v cat
# with PARRALEL
/usr/bin/time R --no-restore --no-save --args FALSE TRUE load $M $N TRUE FALSE <  cghseg.examples.args.R | grep Integration | grep -v cat
# with  PARRALEL+OPTIMIZATION
/usr/bin/time R --no-restore --no-save --args FALSE TRUE load $M $N TRUE TRUE <  cghseg.examples.args.R | grep Integration | grep -v cat

