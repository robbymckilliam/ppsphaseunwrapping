#export JAVA_OPTS="-Xprof -server -Xms1g -Xmx1g"
export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"
scala -nocompdaemon -cp PubSim.jar:Jama-1.0.2.jar:flanagan.jar:colt.jar:RngPack.jar simgaussian.scala
scala -nocompdaemon -cp PubSim.jar:Jama-1.0.2.jar:flanagan.jar:colt.jar:RngPack.jar simwrappedu.scala
