package org.nmdp.ngs.kharsh

/*
 * Sample from a user-supplied distribution.
 *
 * source: http://codereview.stackexchange.com/questions/55409/sampling-from-a-distribution-with-given-probabilities
 * 
 * @author http://codereview.stackexchange.com/users/9607/toto2
 * @version $Id: ProbabilityScript.groovy 22961 2015-06-06 01:12:40Z droe $
 */
public class ProbabilitySampler {
    private Random random = new Random();
    private List<Double> probabilities;
    private double[] cdf;

    public ProbabilitySampler(List<Double> probabilities) {
        super();
        // TODO check sum of probabilities is very close to 1.0
        this.probabilities = Collections.unmodifiableList(probabilities);
        this.cdf = buildCDF(this.probabilities);
    }

    public static double[] buildCDF(List<Double> probabilities) {
        double[] cdf = new double[probabilities.size()];
        cdf[0] = probabilities.get(0);
        for (int i = 1; i < probabilities.size(); i++)
            cdf[i] = cdf[i - 1] + probabilities.get(i);
        return cdf;
    }

    public Integer sample() {
        int index = Arrays.binarySearch(cdf, random.nextDouble());
        return (index >= 0) ? index : (-index - 1);
    }


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
}

