import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;

public class GeneticAlgorithmConfiguration {
    private Random rand;
    private int iterations;
    private double mutation;
	private double initMutation;
	private double crossover;
	private int convergence;
	private double similarity;
    private ChromosomeFitnessValues chromosomePopulation;
    private int numAgents;
    private int size;
    private ChromosomeTimePair bestChromosomeOverall;
    private boolean allowDuplicateGenes;
    // key = chromosome, value = time
	private Map<String, Double> bestChromosomesMap;
	// key = chromosome, value = iteration
	private Map<String, Integer> iterationChromosomeMap;
    private boolean batchProcessing;
    private String workingPath;
    private static final String FITNESS_FILE_NAME = "fitness.txt";

    public GeneticAlgorithmConfiguration() {
        bestChromosomesMap = new HashMap<>();
        iterationChromosomeMap = new HashMap<>();
        chromosomePopulation = new ChromosomeFitnessValues();
        allowDuplicateGenes = false;
        batchProcessing = false;
        rand = new Random();
    }

    public File getFitnessFile() {
        return new File(workingPath + "\\" + FITNESS_FILE_NAME);
    }

    public boolean duplicateGenesEnabled() {
        return allowDuplicateGenes;
    }

    public boolean batchModeEnabled() {
        return batchProcessing;
    }

	public Iterable<String> getBestChromosomes() {
		return bestChromosomesMap.keySet();
	}
	
	public Double getBestChromosomeTime(String chromosome) {
		return bestChromosomesMap.get(chromosome);
	}
	
	public Map<String, Double> getBest() {
		return bestChromosomesMap;
	}
	
	public Map<String, Integer> getChromosomeIterations() {
		return iterationChromosomeMap;
	}
	
	public Integer getChromosomeIteration(String chromosome) {
		return iterationChromosomeMap.get(chromosome);
	}
	
	/**
	 * @return the iterations
	 */
	public int getIterations() {
		return iterations;
	}

	/**
	 * @param iterations the iterations to set
	 */
	public void setIterations(int iterations) {
		this.iterations = iterations;
	}

	/**
	 * @return the mutation
	 */
	public double getMutation() {
		return mutation;
	}

	/**
	 * @param mutation the mutation to set
	 */
	public void setMutation(double mutation) {
		this.mutation = mutation;
	}

	/**
	 * @return the crossover
	 */
	public double getCrossover() {
		return crossover;
	}

	/**
	 * @param crossover the crossover to set
	 */
	public void setCrossover(double crossover) {
		this.crossover = crossover;
	}

	/**
	 * @return the convergence
	 */
	public int getConvergence() {
		return convergence;
	}

	/**
	 * @param convergence the convergence to set
	 */
	public void setConvergence(int convergence) {
		this.convergence = convergence;
	}

	/**
	 * @return the similarity
	 */
	public double getSimilarity() {
		return similarity;
	}

	/**
	 * @param similarity the similarity to set
	 */
	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}

	/**
	 * @return the numAgents
	 */
	public int getNumAgents() {
		return numAgents;
	}

	/**
	 * @param numAgents the numAgents to set
	 */
	public void setNumAgents(int numAgents) {
		this.numAgents = numAgents;
	}

	/**
	 * @return the size
	 */
	public int getSize() {
		return size;
	}

	/**
	 * @param size the size to set
	 */
	public void setSize(int size) {
		this.size = size;
	}

	/**
	 * @return the initMutation
	 */
	public double getInitMutation() {
		return initMutation;
	}

	/**
	 * @param initMutation the initMutation to set
	 */
	public void setInitMutation(double initMutation) {
		this.initMutation = initMutation;
		this.mutation = initMutation;
	}

	public String getBestOverallChromosome() {
        if (bestChromosomeOverall == null)
            return "";

		return bestChromosomeOverall.getChromosome();
	}

	public void addBestChromosome(String chromosome, Double bestTime) {
		// do not overwrite the original time seen for this chromosome
		if(bestChromosomesMap.get(chromosome) == null) {
			this.bestChromosomesMap.put(chromosome, bestTime);

            updateBestOverall(chromosome, bestTime);
        }
	}

    public void updateBestOverall(String chromosome, Double bestTime) {
        if(bestChromosomeOverall == null || bestChromosomeOverall.getTime() > bestTime)
			bestChromosomeOverall = new ChromosomeTimePair(chromosome, bestTime);
	}

    public void increaseMutationRate(Double v) {
        this.mutation += v;
    }

    public void allowDuplicates() {
        allowDuplicateGenes = true;
    }

    public void batchMode() {
        batchProcessing = true;
    }

    public void setWorkingPath() {
        Scanner scan = new Scanner(System.in);
        System.out.println("Provide the path where the simulation is expecting input and provides output:");
        workingPath = scan.nextLine();
        scan.close();
    }

    public String getWorkingPath() {
        return workingPath;
    }

    public void setPopulationFitnessValues() {
        chromosomePopulation.updateFitnessValues();
    }

    public Integer getWinner() {
        return chromosomePopulation.getBestChromosome();
    }

    public Double getWinningFitnessValue() {
        return chromosomePopulation.getWinningFitnessValue();
    }

    public Integer getLoser() {
        return chromosomePopulation.getLoser();
    }

    public void doCrossover() {
        for (int i = 0; i < getSize(); i++) {
            if (rand.nextDouble() < getCrossover()) {
                reproduction(getWinners(), getLoser(), allowDuplicates());
            }
        }
    }

    private class ChromosomeTimePair {
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

    // if this winner has not been seen, add it to the bestChromosomeMap
    // -- returns true if the map has been updated, false otherwise
    boolean updateWinner(int loop) {
        boolean updated = false;
        String chromString = chromosomePopulation.getWinnerString();
        if (bestChromosomesMap.get(chromString) == null) {
            System.out.println("adding to map : " + chromString + "(" + chromosomePopulation.getWinningFitnessValue() + ") : loop = " + loop);
            bestChromosomesMap.put(chromString, chromosomePopulation.getWinningFitnessValue());
            iterationChromosomeMap.put(chromString, (loop + 1));
            updated = true;
        }

        return updated;
    }

    private static int[] getWinners(Double winnerSimilarityRatio, int[][] population, int convergenceFactor) {
        int numTournament = 2;
        // winning gene positions, the start of binary selection
        int[] winners = new int[numTournament];

        int tempConverge = convergenceFactor;
        double diffCount = 0;
        double tempSimilarity = winnerSimilarityRatio;
        while (diffCount / population[0].length <= tempSimilarity) {
            if (diffCount > 0) {
                tempSimilarity -= 0.01;
                tempConverge--;
            }
            diffCount = 0.0;

            int lastWinner = -1;
            for (int i = 0; i < numTournament; i++) {
                winners[i] = selectParentUsingBinaryTournament(lastWinner);
                lastWinner = winners[i];
            }

            // check similarity, if 2 parents are near identical, then the next
            // child will probably duplicate its parent
            // TL;DR - lets prevent in-breeding
            int[] winner1 = population[winners[0]];
            int[] winner2 = population[winners[1]];
            for (int i = 0; i < winner1.length; i++) {
                if (winner1[i] != winner2[i])
                    diffCount++;
            }
        }

        return winners;
    }

    public static int selectParentUsingBinaryTournament(int lastWinner, Integer popSize, int[] popmember) {
        int x = randomNumber(0, popSize);
        int y = randomNumber(0, popSize);
        while (x == y || x == lastWinner || y == lastWinner) {
            x = randomNumber(0, popSize);
            y = randomNumber(0, popSize);
        }

        if (popmember[x] < popmember[y]) {
            return x;
        } else {
            return y;
        }
    }

    public static int randomNumber(int min, int max) {
        double d = min + rand.nextDouble() * (max - min);
        return (int) d;
    }

    private class ChromosomeFitnessValues {
        private double[] fitnessValues;
        private int[][] population;
        private int bestChromosome;
        private int worstChromosome;

        ChromosomeFitnessValues() {
            bestChromosome = -1;
            worstChromosome = -1;
        }

        public String getWinnerString() {
            StringBuilder winnerString = new StringBuilder();

            for(int i = 0; i < population.length ; i++) {
                winnerString.append(population[i]);

                if(i != population.length - 1)
                    winnerString.append(",");
            }

            return winnerString.toString();
        }

        public double[] getFitnessValues() {
            return fitnessValues;
        }

        public int[][] getPopulation() {
            return population;
        }

        public void updateFitnessValues() {
            double largestFitness = Double.MIN_VALUE;
            double prevLarge = largestFitness;

            double smallestFitness = Double.MAX_VALUE;
            double prevSmall = smallestFitness;
            try (Scanner fitnessScanner = new Scanner(getFitnessFile())) {
                for (int i = 0; i < getSize(); i++) {
                    fitnessValues[i] = fitnessScanner.nextDouble();

                    largestFitness = updateWorst(largestFitness, fitnessValues[i]);
                    if(prevLarge != largestFitness) worstChromosome = i;
                    prevLarge = largestFitness;

                    smallestFitness = updateBest(smallestFitness, fitnessValues[i]);
                    if(prevSmall != smallestFitness) bestChromosome = i;
                    prevSmall = smallestFitness;
                }
            } catch(FileNotFoundException e) {
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

        public Integer getBestChromosome() {
            return bestChromosome;
        }

        public Double getWinningFitnessValue() {
            return fitnessValues[bestChromosome];
        }

        public Integer getLoser() {
            return worstChromosome;
        }
    }
}
