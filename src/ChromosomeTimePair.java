/**
 * Created by Erik on 9/4/2016.
 */
public class ChromosomeTimePair {
    private String chromosome;
    private Double time;

    ChromosomeTimePair(String s, Double d) {
        chromosome = s;
        time = d;
    }

    String getChromosome() {
        return chromosome;
    }

    Double getTime() {
        return time;
    }
}
