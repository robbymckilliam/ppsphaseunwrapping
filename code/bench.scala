/**
 * Run simulations for the sparse noisy period estimation problem.
 * Various estimators are available.
 */
import pubsim.poly.PolynomialPhaseSignal
import pubsim.poly.HAF
import pubsim.poly.Babai
import pubsim.poly.Mbest
import pubsim.poly.MaximumLikelihood
import pubsim.poly.SphereDecoder
import pubsim.poly.PolynomialPhaseEstimatorInterface
import pubsim.poly.CPF
import pubsim.poly.CPFHAF
import pubsim.poly.ZW
import pubsim.poly.bounds.GaussianCRB
import pubsim.distributions.GaussianNoise
import pubsim.lattices.reduction.None
import pubsim.lattices.reduction.HKZ
import pubsim.lattices.reduction.LLL

val MIN_BENCH_DURATION : Long = 40000000000L; // (20 secs)  
val Ns = List(30, 60, 100, 200, 350, 500);

def npow(x : Int, t : Int) : Double = if(t<=0) 1.0 else scala.math.pow(x,t)
val p3 = (0 to 3).map( k => 0.25/pubsim.Util.factorial(k) ).toArray //3rd order paramaters
val p5 = (0 to 5).map( k => 0.25/pubsim.Util.factorial(k) ).toArray //5th order paramaters

runbench(Ns, sp3, N, iters, () => new HAF(3,N), "benchHAF3")
runbench(Ns, sp3, N, iters, () => new Mbest(3,N), "benchMbest3")
runbench(Ns, sp5, N, iters, () => new HAF(5,N), "benchHAF5")
runbench(Ns, sp5, N, iters, () => new Mbest(5,N), "benchMbest5")

def runbench(N : Seq[Int], params : Array[Double], snrdB : Double, benchtime : Double, estf : () => PolynomialPhaseEstimatorInterface, name : String) {

val snr = scala.math.pow(10.0, snrdB/10.0)
val noise =  new GaussianNoise(0,1.0/snr/2.0) //variance for real and imaginary parts (divide by 2)

val timelist = Ns.map{ N =>

  val siggen =  new PolynomialPhaseSignal(N) //construct a signal generator
siggen.setNoiseGenerator(noise)
siggen.generateRandomParameters(m)
siggen.generateReceivedSignal

val est = estf()

print(" ... Warming up ... ")
var numiters = 0
val warmupstarttime = System.nanoTime
while(System.nanoTime - warmupstarttime < MIN_BENCH_DURATION/2){
   est.estimate(siggen.getReal, siggen.getImag)
   numiters = numiters+2
}
    
print("Benchmarking ... ")
val benchstarttime = System.nanoTime
for( i <- 1 to numiters) est.estimate(siggen.getReal, siggen.getImag)
val tNano = (System.nanoTime - benchstarttime +(numiters)/2) / numiters //copied from Alan Eliasen java BigInteger benchmarks
val tSec = tNano/1000000000.0 //time in seconds per iteration
println(tSec + " seconds per iteration")

//return time in seconds
tSec

}.toList

val file = java.io.FileWriter("data/" + name)
(Ns, timelist).zipped.foreach { (N, secs) =>
  file.write(N.toString.replace('E', 'e') + "\t" + secs.toString.replace('E', 'e') + "\n")
}
file.close

}


