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
import pubsim.poly.CPFHAF
import pubsim.poly.ZW
import pubsim.poly.bounds.AngularLeastSquaresVariance
import pubsim.poly.bounds.GaussianCRB
import pubsim.distributions.GaussianNoise
import pubsim.distributions.circular.ProjectedNormalDistribution
import pubsim.lattices.reduction.None
import pubsim.lattices.reduction.HKZ
import pubsim.lattices.reduction.LLL

val iters = 100
val N = 199
val sp3 = (0 to 3).map( k => 0.5/pubsim.Util.factorial(k)/scala.math.pow(N,k-1) ).toArray //3rd order paramater that work for the HAF and CPF
val bp3 = (0 to 3).map( k => 0.5/pubsim.Util.factorial(k)/2.0 ).toArray //3rd order paramater that do not work for the HAF and CPF
val sp5 = (0 to 5).map( k => 0.5/pubsim.Util.factorial(k)/scala.math.pow(N,k-1) ).toArray //5rd order paramater that work for the HAF
val bp5 = (0 to 5).map( k => 0.5/pubsim.Util.factorial(k)/2.0 ).toArray //5rd order paramater that do not work for the HAF
val sp6 = (0 to 6).map( k => 0.5/pubsim.Util.factorial(k)/scala.math.pow(N,k-1) ).toArray //5rd order paramater that work for the HAF

val starttime = (new java.util.Date).getTime
/*
runsim(-5 to 15, sp3, N, iters, () => new HAF(3,N), "HAFm3small")
runsim(-5 to 15, bp3, N, iters, () => new HAF(3,N), "HAFm3big")
runsim(-5 to 15, sp3, N, iters, () => new CPF(N), "CPFm3small")
runsim(-5 to 15, bp3, N, iters, () => new CPF(N), "CPFm3big")
runsim(-5 to 15, sp3, N, iters, () => new Babai(3,N, new HKZ()), "Babaim3small")
runsim(-5 to 15, bp3, N, iters, () => new Babai(3,N, new HKZ()), "Babaim3big")
runsim(-5 to 15, sp3, N, iters, () => new Mbest(3,N, 20*N, new HKZ()), "Mbestm3small")
runsim(-5 to 15, bp3, N, iters, () => new Mbest(3,N, 20*N, new HKZ()), "Mbestm3big")
runcrb(-5 to 15, 3, N, "crbm3")
runlsuclt(-5 to 15, 3, N, "lsucltm3")

runsim(-5 to 35 by 2, bp3, N, iters, () => new ZW(3,N,N/3,N/3+1), "ZWm3big")
runsim(-5 to 35 by 2, bp3, N, iters, () => new Babai(3,N, new HKZ()), "Babaim3bigrange")
runsim(-5 to 35 by 2, bp3, N, iters, () => new Mbest(3,N, 20*N, new HKZ()), "Mbestm3bigrange")
runcrb(-5 to 35, 3, N, "crbm3range")
runlsuclt(-5 to 35, 3, N, "lsucltm3range")

runsim(5 to 30, sp5, N, iters, () => new HAF(5,N), "HAFm5small")
runsim(5 to 30, sp5, N, iters, () => new CPFHAF(5,N), "CPFHAFm5small")
runsim(5 to 30, sp5, N, iters, () => new Babai(5,N, new LLL()), "Babaim5small")
runsim(5 to 30, sp5, N, iters, () => new Mbest(5,N, 20*N, new LLL()), "Mbestm5small")
runcrb(5 to 30, 5, N, "crbm5")
runlsuclt(5 to 30, 5, N, "lsucltm5")
*/

runsim(5 to 30, sp6, N, iters, () => new HAF(6,N), "HAFm6small")
runsim(5 to 30, sp6, N, iters, () => new CPFHAF(6,N), "CPFHAFm6small")
runsim(5 to 30, sp6, N, iters, () => new Babai(6,N, new LLL()), "Babaim6small")
runsim(5 to 30, sp6, N, iters, () => new Mbest(6,N, 20*N, new LLL()), "Mbestm6small")
runcrb(5 to 30, 6, N, "crbm6")
runlsuclt(5 to 30, 6, N, "lsucltm6")

val runtime = (new java.util.Date).getTime - starttime
println("Simulation finshed in " + (runtime/1000.0) + " seconds.\n")




/** Runs a simulation with given parameters and stores output in a file */
def runsim(snrdbs : Seq[Int], params : Array[Double], N : Int, iters : Int, estf : () => PolynomialPhaseEstimatorInterface, name : String) {

  val m = params.size - 1
  val SNRs = snrdbs.map(db => scala.math.pow(10.0, db/10.0))
  val noises =  SNRs.map( snr => new GaussianNoise(0,1.0/snr/2.0) ) //variance for real and imaginary parts (divide by 2)
  
  print("Running " + name + " ")
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
  val files = (0 to m).map( v => new java.io.FileWriter("data/" + name + "p" + v) ) //list of files to write to
  (mselist, snrdbs).zipped.foreach{ (mse, snr) =>
    for ( i <- files.indices ) 
      files(i).write(snr.toString.replace('E', 'e') + "\t" + (mse(i)/iters).toString.replace('E', 'e')  + "\n") 
  }
  for( f <- files ) f.close //close all the files we wrote to 
  
  val estruntime = (new java.util.Date).getTime - eststarttime
  println(" finished in " + (estruntime/1000.0) + " seconds.")

}


//Gaussian CRB plots
def runcrb(snrdbs : Seq[Int], m : Int, N : Int, name : String) {
  val SNRs = snrdbs.map(db => scala.math.pow(10.0, db/10.0))
  val noises =  SNRs.map( snr => new GaussianNoise(0,1.0/snr/2.0) ) //variance for real and imaginary parts (divide by 2)
    print("Computing " + name + " ")
    val gausscrb = new GaussianCRB(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("data/" + name + "p" + v) )
    for ( i <- files.indices ) {
      for( s <- snrdbs.indices ) {
	val v = noises(s).getVariance
	files(i).write(snrdbs(s).toString.replace('E', 'e') + "\t" + gausscrb.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
      print(".")
    }
    for ( f <- files ) f.close //close all the files we wrote to 
    println(" done")
  }


//Least squares unwrapping asymptotics
def runlsuclt(snrdbs : Seq[Int], m : Int, N : Int, name : String) {
  val SNRs = snrdbs.map(db => scala.math.pow(10.0, db/10.0))
  val noises =  SNRs.map( snr => new GaussianNoise(0,1.0/snr/2.0) ) //variance for real and imaginary parts (divide by 2)
  print("Computing " + name + " ")
    val lsuasymp = new AngularLeastSquaresVariance(N,m)
    val files = (0 to m).map( v => new java.io.FileWriter("data/" + name + "p" + v) )
    for ( i <- files.indices ) {
      for( s <- snrdbs.indices ) {
	val v = new ProjectedNormalDistribution(0,noises(s).getVariance)
	files(i).write(snrdbs(s).toString.replace('E', 'e') + "\t" + lsuasymp.getBound(i,v).toString.replace('E', 'e') + "\n")
      }
      print(".")
    }
    for ( f <- files ) f.close //close all the files we wrote to 
    println(" done")
  }


