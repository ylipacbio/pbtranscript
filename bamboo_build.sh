#!/bin/bash -e
rm -rf prebuilt
rm -rf deployment
mkdir -p prebuilt
PBBAM=`/bin/ls -t tarballs/pbbam*-x86_64.tgz|head -1`
HTSLIB=`/bin/ls -t tarballs/htslib-*tgz|head -1`
BLASR=`/bin/ls -t tarballs/blasr-*tgz|head -1`
BLASR_LIBCPP=`/bin/ls -t tarballs/blasr_libcpp*tgz|head -1`
PBDAGCON=`/bin/ls -t tarballs/pbdagcon-*tgz|head -1`

HAVE_NCURSES=$PWD/prebuilts/`basename $NCURSES .tgz`
HAVE_HTSLIB=$PWD/prebuilts/`basename $HTSLIB .tgz`
HAVE_BLASR=$PWD/prebuilts/`basename $BLASR .tgz`
HAVE_BLASR_LIBCPP=$PWD/prebuilts/`basename $BLASR_LIBCPP .tgz`
HAVE_PBBAM=$PWD/prebuilts/`basename $PBBAM -x86_64.tgz`
HAVE_PBDAGCON=$PWD/prebuilts/`basename $PBDAGCON .tgz`
mkdir -p  \
         $HAVE_HTSLIB \
         $HAVE_BLASR_LIBCPP \
         $HAVE_BLASR \
         $HAVE_PBBAM \
         $HAVE_PBDAGCON
tar zxf $HTSLIB -C $HAVE_HTSLIB
tar zxf $PBBAM -C prebuilts
tar zxf $BLASR -C $HAVE_BLASR
tar zxf $BLASR_LIBCPP -C $HAVE_BLASR_LIBCPP
tar zxf $PBDAGCON -C $HAVE_PBDAGCON

type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load git/2.8.3
module load gcc/4.9.2
module load ccache/3.2.3
module load graphviz

cat > pitchfork/settings.mk << EOF
DISTFILES             = ${PWD}/.distfiles
CCACHE_DIR            = ${PWD}/.ccache
PIP_CACHE             = ${PWD}/.pip
# from Herb
HAVE_OPENSSL      = /mnt/software/o/openssl/1.0.2a
HAVE_PYTHON       = /mnt/software/p/python/2.7.9/bin/python
HAVE_BOOST        = /mnt/software/b/boost/1.58.0
HAVE_ZLIB         = /mnt/software/z/zlib/1.2.8
HAVE_SAMTOOLS     = /mnt/software/s/samtools/1.3.1mobs
HAVE_NCURSES      = /mnt/software/n/ncurses/5.9
# from MJ
HAVE_HDF5         = /mnt/software/a/anaconda2/4.2.0
HAVE_OPENBLAS     = /mnt/software/o/openblas/0.2.14
HAVE_CMAKE        = /mnt/software/c/cmake/3.2.2/bin/cmake
#
HAVE_HTSLIB           = $HAVE_HTSLIB
HAVE_BLASR_LIBCPP     = $HAVE_BLASR_LIBCPP
HAVE_BLASR            = $HAVE_BLASR
HAVE_PBDAGCON         = $HAVE_PBDAGCON

pbtranscript_REPO     = $PWD/repos/pbtranscript

pbbam_REPO            = ssh://git@bitbucket.nanofluidics.com:7999/sat/pbbam.git
ConsensusCore_REPO    = ssh://git@bitbucket.nanofluidics.com:7999/sat/ConsensusCore.git
bam2fastx_REPO        = ssh://git@bitbucket.nanofluidics.com:7999/sat/bam2fastx.git
pbcoretools_REPO      = ssh://git@bitbucket.nanofluidics.com:7999/sat/pbcoretools.git
pbdagcon_REPO         = ssh://git@bitbucket.nanofluidics.com:7999/sat/pbdagcon.git
daligner_REPO         = ssh://git@bitbucket.nanofluidics.com:7999/sat/daligner.git
dazzdb_REPO           = ssh://git@bitbucket.nanofluidics.com:7999/sat/dazz_db.git
GenomicConsensus_REPO = ssh://git@bitbucket.nanofluidics.com:7999/sat/GenomicConsensus.git
pbcopper_REPO         = ssh://git@bitbucket.nanofluidics.com:7999/sat/pbcopper.git
pbcore_REPO           = ssh://git@bitbucket.nanofluidics.com:7999/sat/pbcore.git
pblaa_REPO            = ssh://git@bitbucket.nanofluidics.com:7999/sat/pblaa.git
seqan_REPO            = ssh://git@bitbucket.nanofluidics.com:7999/sat/seqan.git
unanimity_REPO        = ssh://git@bitbucket.nanofluidics.com:7999/sat/unanimity.git

pbcommand_REPO        = ssh://git@bitbucket.nanofluidics.com:7999/sl/pbcommand.git
EOF
# this doesn't work because of pysam
# make -C pitchfork -j15 pbtranscript
cd pitchfork
make -j15 pbtranscript
# disable SGE
ln -sfn /bin/false deployment/bin/qstat
