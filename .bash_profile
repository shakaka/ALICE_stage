
if [ -f "$HOME/.bashrc" ] ; then
    source $HOME/.bashrc
fi

#ALIROOT define in $WORKDIR/code/local/script/alice-alienx-env.sh
#version_root=v5-28-00a
#version_geant=v1-11
#version_aliroot=v4-21-18-AN

#LD_LIBRARY_PATH=/usr/lib:/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/cynthia

#WORK ON IPNGRILLE
if [ -f "/ipn/storage1/env/gridenv.sh" ] ; then
    source /ipn/storage1/env/gridenv.sh
    export PATH=${PATH}:/nfs1/ipno_soft/cmake/2.8.4/bin

    export INSTALLDIR=/ipn/storage1/alice/users/cynthia

    export WORKDIR=/nfs1/alice/users/cynthia
    export SCRATCH=/nfs1/scratch/cynthia
    #export LD_LIBRARY_PATH=$WORKDIR/code/CGC/CGC_Gelis_Code_ompi/gsl_x86-64/lib/:$LD_LIBRARY_PATH

    #source /nfs1/env/gridenv_pgi7.1.sh
    #export PATH=/nfs1/openmpi/1.2.6/pgi_x86_64_without-tm/bin:${PATH}
    #export LD_LIBRARY_PATH=/nfs1/openmpi/1.2.6/pgi_x86_64_without-tm/lib:${LD_LIBRARY_PATH}

    hard_platform=`uname  -i`
#    if [ "${hard_platform}" = "x86_64" ]; then
#	export LD_LIBRARY_PATH=${PGI}/linux86-64/7.1/libso:${LD_LIBRARY_PATH}
#    else
#	export LD_LIBRARY_PATH=${PGI}/linux86/7.1/libso:${LD_LIBRARY_PATH}
#    fi
else 
#WORK Environment
    #export WORKDIR=/projet/alice/Cynthia/
    #export SCRATCH=/projet/alice/Cynthia/temp
    #for ipnalice9
    export WORKDIR=/home/cynthia/alice
    export SCRATCH=/vol0/cynthia/scratch
fi




#CASCADE Environment
export CERN_LIBS="/import/cern/99/lib -lmathlib -lkernlib -lpacklib"
export PYTHIA_LIBS="/import/cern/99/lib -l pythia6125"
export HZTOOL_LIBS="/nfs1/alice/users/cynthia/code/hztool/lib -lhztool"
#ALICE Environment

#PATH
export PATH=${PATH}:/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11R6/bin

#MAN
export MANPATH=${MANPATH}:/usr/share/man

return


#ROOT 
if [ -f "/ipn/storage1/env/gridenv.sh" ] ; then
    export ALICE=/nfs1/alice/users/cynthia/code/
    export ROOTSYS=$ALICE/root/${version_root}/
    export PATH=${PATH}:${ROOTSYS}/bin
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
    export MANPATH=$MANPATH:$ROOTSYS/man
fi

# Skip the rest if ROOT hasn't been compiled yet
which root-config > /dev/null 2> /dev/null
if [ $? != 0 ]; then
  echo ""
  echo " * ROOT has not been installed yet. Re-source the needed environment variables"
  echo "   after building ROOT!"
  echo ""
  return
fi

#export ALICE_TARGET=`root-config --arch`
#G3
G3LIB=${ALICE}/geant3/${version_geant}/lib/tgt_${ALICE_TARGET}
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G3LIB}/

#AliRoot
ALICE_LEVEL=aliroot/${version_aliroot}
export ALICE_ROOT=${ALICE}/${ALICE_LEVEL}
export ALICE_ARCH=`uname`
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}
PATH=${PATH}:${ALICE_ROOT}/bin/tgt_${ALICE_TARGET}

#Alien
#GSHELL_GCC=`which gcc`
#ALIEN=$ALICE/alien
#GSHELL_ROOT=${ALIEN}/api/bin
#GLOBUS_LOCATION=${ALIEN}/globus
#following line should be commented if using GRIF
#X509_CERT_DIR=${ALIEN}/globus/share/certificates
#PATH=${GSHELL_ROOT}:${GLOBUS_LOCATION}/bin:${PATH}
#export LD_LIBRARY_PATH=${ALIEN}/api/lib:${GLOBUS_LOCATION}/lib:${LD_LIBRARY_PATH}

export PATH=$PATH:.:~$LOGNAME/bin

