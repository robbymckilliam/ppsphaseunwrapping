/**
 * Run simulations for the sparse noisy period estimation problem.
 * Various estimators are available.
 */
import pubsim.poly.PolynomialPhaseSignal
import pubsim.poly.KitchenEstimator
import pubsim.poly.DPTEstimator
import pubsim.poly.MbestEstimator
import pubsim.poly.BabaiEstimator
import pubsim.poly.MbestEstimator
import pubsim.poly.MaximumLikelihood
import pubsim.poly.SphereDecoderEstimator
import pubsim.poly.PolynomialPhaseEstimator
import pubsim.poly.CubicPhaseFunction
import pubsim.poly.bounds.AngularLeastSquaresVariance
import pubsim.poly.bounds.GaussianCRB
import pubsim.distributions.GaussianNoise
import pubsim.distributions.circular.ProjectedNormalDistribution

val iters = 100 //number of Monte-Carlo trials.
val Ns = List(39)//,249) //values of N we will generate curves for
val ms = List(3) //order of our polynomial phase signals

def gp( a : Double, r : Double, n : Int ) : Double = a * scala.math.pow(r,n)

//returns an array of noise distributions with a logarithmic scale
val SNRdBs = 0 to 21
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
def noises =  SNRs.map( snr => new GaussianNoise(0,1.0/snr/2.0) ) //variance for real and imaginary parts (divide by 2)

//Returns a list of functions that return estimators we will run (factory patern to enable parallelism)
def estfactory(m : Int, N : Int) : List[() => PolynomialPhaseEstimator] = {
  var ret = List( 
    () => new KitchenEstimator(m,N),
    () => new DPTEstimator(m,N),
    //() => new BabaiEstimator(m,N),
    () => new MbestEstimator(m,N,4*N) 
  )
  //add the sphere decoder and Least squares estimators if N and m are small
  if( N < 60 ) ret = ret :+ ( () => new SphereDecoderEstimator(m,N) )
  //if( N < 60 && m <= 2 ) ret = ret :+ ( () => new MaximumLikelihood(m,N) )
if( m==3 ) ret = ret :+ ( () => new CubicPhaseFunction(N) )
  return ret
}

val starttime = (new java.util.Date).getTime

//for all the the values of N and m.
for( N <-  Ns; m <- ms ) {

  for(estf <- estfactory(m,N) ){
    
    val estname = estf().getClass.getSimpleName + "N" + N.toString + "m" + m + "Gaussian"
    print("Running " + estname)
    val eststarttime = (new java.util.Date).getTime
    
    //for all the noise distributions (in parallel threads)
    val mselist = noises.par.map { noise =>

      val siggen =  new PolynomialPhaseSignal(N) //construct a signal generator
      siggen.setNoiseGenerator(noise)
      val est = estf() //construct an estimator
				  
      var mse = new Array[Double](m+1) //storage for the mses
      for( itr <- 1 to iters ) {
	siggen.generateRandomParameters(m)
	val p0 = siggen.getParameters
	siggen.generateReceivedSignal
	val err = est.error(siggen.getReal, siggen.getImag, p0)
	for( i <- mse.indices ) mse(i) += err(i)*err(i)
      }
      print(".")
      mse //last thing is what gets returned 
    }.toList
    
    //now write all the data to a file
    val files = (0 to m).map( v => new java.io.FileWriter(estname + "p" + v) ) //list of files to write to
    (mselist, noises).zipped.foreach{ (mse, noise) =>
      for ( i <- files.indices ) 
	files(i).write(noise.getVariance.toString.replace('E', 'e') + "\t" + (mse(i)/iters).toString.replace('E', 'e')  + "\n") 
    }
    for( f <- files ) f.close //close all the files we wrote to 
    
    val estruntime = (new java.util.Date).getTime - eststarttime
    println(" finished in " + (estruntime/1000.0) + " seconds.")

  }
  
    //Gaussian CRB plots
  {
    val gausscrb = new GaussianCRB(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("GaussCRBN" + N + "m" + m  + "p" + v) )
    for ( i <- files.indices ) {
      for(noise <- noises ) {
	val v = noise.getVariance
	files(i).write(v.toString.replace('E', 'e') + "\t" + gausscrb.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
    }
    for ( f <- files ) f.close //close all the files we wrote to 
  }

  //Least squares unwrapping asymptotics
  {
    val lsuasymp = new AngularLeastSquaresVariance(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("LSUasympGaussianN" + N + "m" + m  + "p" + v) )
    for ( i <- files.indices ) {
      for(noise <- noises ) {
	val v = new ProjectedNormalDistribution(0,noise.getVariance)
	files(i).write(noise.getVariance.toString.replace('E', 'e') + "\t" + lsuasymp.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
    }
    for ( f <- files ) f.close //close all the files we wrote to 
  }

}

val runtime = (new java.util.Date).getTime - starttime
println("Simulation finshed in " + (runtime/1000.0) + " seconds.\n")
