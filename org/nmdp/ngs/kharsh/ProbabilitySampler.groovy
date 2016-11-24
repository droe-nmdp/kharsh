package org.nmdp.ngs.kharsh

/*
 * Sample from a user-supplied distribution.
 *
 * source: http://codereview.stackexchange.com/questions/55409/sampling-from-a-distribution-with-given-probabilities
 * 
 * @author http://codereview.stackexchange.com/users/9607/toto2
 * @author Dave Roe
 * @version $Id: ProbabilitySampler.groovy 23632 2015-07-04 00:34:47Z droe $
 */
public class ProbabilitySampler {
    private Random random = new Random();
    private List<Double> probabilities;
    private double[] cdf;

    public ProbabilitySampler(List<Double> p) {
        super();

        probabilities = p.collect()
        // normalize
        double total = 0;
        for (int i = 0; i < probabilities.size(); i++) { 
            total += probabilities.get(i)
        }
        for (int i = 0; i < probabilities.size(); i++) { 
            probabilities.set(i, new Double(probabilities.get(i) / total))
        }
        
        this.probabilities = Collections.unmodifiableList(probabilities);
        this.cdf = buildCDF(this.probabilities);
    }

    public static double[] buildCDF(List<Double> probabilities) {
        double[] cdf = new double[probabilities.size()];        
        cdf[0] = probabilities.get(0);
        for (int i = 1; i < probabilities.size(); i++) { 
            cdf[i] = cdf[i - 1] + probabilities.get(i);
        }
        //System.err.println "buildCDF: cdf[${probabilities.size()-1}]=" + cdf[probabilities.size()-1]//todo: remove(should be 1)
            
        return cdf;
    }

    public Integer sample() {
        int index = Arrays.binarySearch(cdf, random.nextDouble());
        Integer ret = (index >= 0) ? index : (-index - 1);
        //System.err.println "cdf=${cdf}"//todo
        //System.err.println "probability=${probabilities.get(ret)}"
        //System.err.println "ret=${ret}"//todo
        return ret
    }

    public double getProbability(Integer index) {
        return probabilities.get(index)
    } // getProbability

/*
    public static void main(String[] args) {
        List<Double> probabilities = Arrays.asList(0.32, 0.68);
        ProbabilitySampler probabilitySampler = new ProbabilitySampler(probabilities);

        int nSamples = 100000;
        final List<Integer> distribution = new ArrayList<>(Collections.nCopies(probabilities.size(), 0));
        IntStream
            .range(0, nSamples)
            .map(i -> probabilitySampler.sample())
            .forEach(randomItem -> distribution.set(randomItem, distribution.get(randomItem) + 1));
        System.out.println(distribution);
    }
*/
}

