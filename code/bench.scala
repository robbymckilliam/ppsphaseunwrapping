/**
 * Run simulations for the sparse noisy period estimation problem.
 * Various estimators are available.
 */

import pubsim.poly.MbestEstimator
import pubsim.poly.PolynomialPhaseSignal
import pubsim.distributions.GaussianNoise

val m = 8;
val Ns = List(30, 60, 100, 200, 350, 500, 1000);

val MIN_BENCH_DURATION : Long = 40000000000L; // (20 secs)  

val snrdB = 7
val snr = scala.math.pow(10.0, snrdB/10.0)
val noise =  new GaussianNoise(0,1.0/snr/2.0) //variance for real and imaginary parts (divide by 2)

for(N <- Ns){

val siggen =  new PolynomialPhaseSignal(N) //construct a signal generator
siggen.setNoiseGenerator(noise)
siggen.generateRandomParameters(m)
siggen.generateReceivedSignal

val est = new MbestEstimator(m,N,4*N) //contruct an estimator to benchmark

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
val tSec = tNano/1000000000.0 //time in milliseconds per iteration
println(tSec + " seconds per iteration")

}
