export JAVA_OPTS="-Xprof -server -Xms1g -Xmx1g"
#export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"

CP=""
for f in java/lib/*.jar
do
CP=$CP:${f}
done

scala -nocompdaemon -cp $CP bench.scala
