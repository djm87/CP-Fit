#! /bin/sh -f
RUN_PATH=$PWD
MAUD_PATH='/mnt/home/thrust2/djm87/soft/Maud/'
cd $MAUD_PATH
MAUD_PATH=`pwd`
./jdk/bin/java -XX:ParallelGCThreads=1 -mx8192M --add-opens java.base/java.net=ALL-UNNAMED -Djava.awt.headless=true -Duser.dir=$MAUD_PATH -cp lib/Maud.jar:lib/ij.jar:lib/jgap.jar:lib/Help.jar:lib/EsquiClient.jar:lib/com.github.tschoonj.xraylib.jar:lib/joone-engine.jar:lib/newt.all.jar:lib/jdic.jar:lib/jdom.jar:lib/sqlite-jdbc.jar:lib/jmol.jar:lib/jgaec.jar:lib/Files.jar:lib/xgridlib.jar:lib/xgridagent.jar:lib/jogl.all.jar:lib/Examples.jar:lib/commons-math.jar:lib/rome.jar:lib/nativewindow.all.jar:lib/Images.jar:lib/swingx.jar:lib/jdic_stub.jar:lib/MySQL-ConnectorJ.jar:lib/HTTPClient.jar:lib/miscLib.jar:lib/gluegen-rt.jar com.radiographema.MaudText -file $RUN_PATH/Maud_Batch_input_"$1".ins
