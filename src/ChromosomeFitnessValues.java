import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * Created by Erik on 9/4/2016.
 *
 * Moved from GeneticAlgorithm configuration into its own class
 */
class ChromosomeFitnessValues {
    private Double[] fitnessValues;
    private Integer[][] population;
    private Integer bestChromosome;
    private Integer worstChromosome;

    ChromosomeFitnessValues() {
        bestChromosome = -1;
        worstChromosome = -1;
    }

    String getWinnerString() {
        StringBuilder winnerString = new StringBuilder();

        for(int i = 0; i < population.length ; i++) {
            winnerString.append(population[i]);

            if(i != population.length - 1)
                winnerString.append(",");
        }

        return winnerString.toString();
    }

    Double[] getFitnessValues() {
        return fitnessValues;
    }

    Integer[][] getPopulation() {
        return population;
    }

    void updateFitnessValues(File fitnessFileName) {
        double largestFitness = Double.MIN_VALUE;
        double prevLarge = largestFitness;

        double smallestFitness = Double.MAX_VALUE;
        double prevSmall = smallestFitness;
        try (Scanner fitnessScanner = new Scanner(fitnessFileName)) {
            for (int i = 0; i < fitnessValues.length; i++) {
                fitnessValues[i] = fitnessScanner.nextDouble();

                largestFitness = updateWorst(largestFitness, fitnessValues[i]);
                if(prevLarge != largestFitness) worstChromosome = i;
                prevLarge = largestFitness;

                smallestFitness = updateBest(smallestFitness, fitnessValues[i]);
                if(prevSmall != smallestFitness) bestChromosome = i;
                prevSmall = smallestFitness;
            }
        } catch(FileNotFoundException e) {
            // do nothing
        }
    }

    private double updateWorst(double largestFitness, double currentFitness) {
        if (currentFitness > largestFitness)
            largestFitness = currentFitness;

        return largestFitness;
    }

    private double updateBest(double smallFitness, double currentFitness) {
        if (currentFitness > smallFitness)
            smallFitness = currentFitness;

        return smallFitness;
    }

    Integer getBestChromosome() {
        return bestChromosome;
    }

    Double getWinningFitnessValue() {
        return fitnessValues[bestChromosome];
    }

    Integer getLoser() {
        return worstChromosome;
    }
}
