#export JAVA_OPTS="-Xprof -server -Xms1g -Xmx1g"
export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"
scala -nocompdaemon -cp java/lib/PubSim.jar:java/lib/Jama-1.0.3.jar:java/lib/flanagan.jar:java/lib/colt.jar:java/lib/RngPack.jar:java/lib/Poly.jar simgaussian.scala
