/**
 * Run simulations of various polynomial phase estimators.
 */
import pubsim.poly.PolynomialPhaseSignal
import pubsim.poly.HAF
import pubsim.poly.Babai
import pubsim.poly.Mbest
import pubsim.poly.MaximumLikelihood
import pubsim.poly.SphereDecoder
import pubsim.poly.PolynomialPhaseEstimatorInterface
import pubsim.poly.CPF
import pubsim.poly.ZW
import pubsim.poly.bounds.AngularLeastSquaresVariance
import pubsim.poly.bounds.GaussianCRB
import pubsim.distributions.GaussianNoise
import pubsim.distributions.circular.ProjectedNormalDistribution
import pubsim.lattices.reduction.None
import pubsim.lattices.reduction.HKZ
import pubsim.lattices.reduction.LLL

val iters = 500 //number of Monte-Carlo trials.
val Ns = List(199) //values of N we will generate curves for
val ms = List(3) //order of our polynomial phase signals

//Returns a list of functions that return estimators we will run (factory patern to enable parallelism)
def estfactory(m : Int, N : Int) : List[() => PolynomialPhaseEstimatorInterface] = {
  var ret = List( 
    //() => new Kitchen(m,N),
    () => new HAF(m,N),
    //() => new Babai(m,N, if(m<=3) new HKZ() else new LLL()),
    //() => new Mbest(m,N,4*N, if(m<=3) new HKZ() else new LLL()),
    () => new ZW(m,N,N/m,N/m+1)
  )
  //add the sphere decoder and Least squares estimators if N and m are small
  if( N < 60 ) ret = ret :+ ( () => new SphereDecoder(m,N) )
  //if( N < 60 && m <= 2 ) ret = ret :+ ( () => new MaximumLikelihood(m,N) )
if( m==3 ) ret = ret :+ ( () => new CPF(N) )
  return ret
}

val starttime = (new java.util.Date).getTime

//for all the the values of N and m.
for( N <-  Ns; m <- ms ) {

//returns an array of noise distributions with a logarithmic scale
val SNRdBs = if(m==3) -5 to 35 else 5 to 30
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
def noises =  SNRs.map( snr => new GaussianNoise(0,1.0/snr/2.0) ) //variance for real and imaginary parts (divide by 2)

  //parameters in the range of DPT/HAF/CPF etc
  def params = (0 to m).map( k => 0.5/pubsim.Util.factorial(k)/scala.math.pow(N,k-1) ).toArray

  for(estf <- estfactory(m,N) ){
    
    val estname = estf().getClass.getSimpleName + "N" + N.toString + "m" + m
    print("Running " + estname + " ")
    val eststarttime = (new java.util.Date).getTime
    
    //for all the noise distributions (in parallel threads)
    val mselist = noises.par.map { noise =>

      val siggen =  new PolynomialPhaseSignal(N) //construct a signal generator
      siggen.setNoiseGenerator(noise)
      siggen.setParameters(params)
      val est = estf() //construct an estimator
				  
      var mse = new Array[Double](m+1) //storage for the mses
      for( itr <- 1 to iters ) {
	//siggen.generateRandomParameters(m)
	val p0 = siggen.getParameters
	siggen.generateReceivedSignal
	val err = est.error(siggen.getReal, siggen.getImag, p0)
	for( i <- mse.indices ) mse(i) += err(i)*err(i)
      }
      print(".")
      mse //last thing is what gets returned 
    }.toList
    
    //now write all the data to a file
    val files = (0 to m).map( v => new java.io.FileWriter("data/" + estname + "p" + v) ) //list of files to write to
    (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
      for ( i <- files.indices ) 
	files(i).write(snr.toString.replace('E', 'e') + "\t" + (mse(i)/iters).toString.replace('E', 'e')  + "\n") 
    }
    for( f <- files ) f.close //close all the files we wrote to 
    
    val estruntime = (new java.util.Date).getTime - eststarttime
    println(" finished in " + (estruntime/1000.0) + " seconds.")

  }
  
    //Gaussian CRB plots
  {
    print("Computing Gaussian CRBS ")
    val gausscrb = new GaussianCRB(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("data/GaussCRBN" + N + "m" + m  + "p" + v) )
    for ( i <- files.indices ) {
      for( s <- SNRdBs.indices ) {
	val v = noises(s).getVariance
	files(i).write(SNRdBs(s).toString.replace('E', 'e') + "\t" + gausscrb.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
      print(".")
    }
    for ( f <- files ) f.close //close all the files we wrote to 
    println(" done")
  }

  //Least squares unwrapping asymptotics
  {
    print("Computing LSU asymptotic variance ")
    val lsuasymp = new AngularLeastSquaresVariance(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("data/LSUasympGaussianN" + N + "m" + m  + "p" + v) )
    for ( i <- files.indices ) {
      for( s <- SNRdBs.indices ) {
	val v = new ProjectedNormalDistribution(0,noises(s).getVariance)
	files(i).write(SNRdBs(s).toString.replace('E', 'e') + "\t" + lsuasymp.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
      print(".")
    }
    for ( f <- files ) f.close //close all the files we wrote to 
    println(" done")
  }

}

val runtime = (new java.util.Date).getTime - starttime
println("Simulation finshed in " + (runtime/1000.0) + " seconds.\n")
