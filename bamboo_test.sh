#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load git/2.8.3
module load gcc/4.9.2
module load ccache/3.2.3
module load graphviz

cat > pbtranscript_dummy.xml << EOF
<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="nosetests" tests="2" errors="0" failures="0" skip="0">
  <testcase classname="dummy.system" name="pwd" time="0.00">
    <system-out><![CDATA[`pwd`]]></system-out>
  </testcase>
  <testcase classname="dummy.system" name="hostname" time="0.00">
    <system-out><![CDATA[`hostname`]]></system-out>
  </testcase>
</testsuite>
EOF
source pitchfork/deployment/setup-env.sh
export PYTHONWARNINGS="ignore"
nosetests --verbose --with-xunit --xunit-file=pbtranscript_nose.xml \
    ../repos/pbtranscript/tests/unit/*.py

chmod +w -R repos/pbtranscript
