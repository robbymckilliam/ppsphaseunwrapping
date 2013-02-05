/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly.bounds;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Vector;

/**
 * Compute the modified CRB for the oversampled DPT estimator
 * @author Robby McKilliam
 */
public class OversampledDPTCRB extends GaussianCRB {

    /** The factor that the DPT must be oversampled */
    private final double dptos;

    public OversampledDPTCRB(int N, int m){
        super(N,m);
        dptos = Math.pow( (N*1.0)/m , (m-1.0)/(m+1.0));
    }

    @Override
    public double getBound(int m, double variance){
        return Math.pow(dptos, 2*m)*variance/
                    (4*Math.PI * Math.PI * Math.pow(N,2*m+1)) * C.get(m, m);
    }

    //** Write CRB data to a file */
    public static void main(String[] args) throws Exception  {

        int N = 256;
        int m = 3;

        OversampledDPTCRB bound = new OversampledDPTCRB(N, m);

        double from_log_snr = 20.0;
        double to_log_snr = -3;
        double step_log_snr = -0.5;

        Vector<Double> snr_array = new Vector<Double>();
        Vector<Double> snr_db_array = new Vector<Double>();
        for(double snrdb = from_log_snr; snrdb >= to_log_snr; snrdb += step_log_snr){
            snr_db_array.add(new Double(snrdb));
            snr_array.add(new Double(Math.pow(10.0, ((snrdb)/10.0))));
        }

        for(int j = 0; j <= m; j++){
            String fname = bound.getClass().getName();
            fname = fname.concat(Integer.toString(N) + "_p" + j);
            File file = new File(fname);
            BufferedWriter writer =  new BufferedWriter(new FileWriter(file));
            for(int i = 0; i < snr_array.size(); i++){
                double var = 0.5/snr_array.get(i);
                writer.write(
                        (new Double(var)).toString().replace('E', 'e')
                        + "\t" + (new Double(bound.getBound(j, var))).toString().replace('E', 'e'));
                writer.newLine();
            }
            writer.close();
        }

    }

    

}
