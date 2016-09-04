import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

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
                reproduction(getWinners(getSimilarity(), chromosomePopulation.getPopulation(), getConvergence()),
                    getLoser(), duplicateGenesEnabled(), chromosomePopulation.getPopulation());
            }
        }
    }

    private void reproduction(int[] winners, int loser, boolean allowDup, int[][] population) {
        int[] child;
        if (allowDup) {
            child = orderOneCrossover(population[winners[0]], population[winners[1]]);
        } else {
            child = orderOneCrossoverWoRepetition(population[winners[0]], population[winners[1]]);
        }
        boolean duplicate = duplicate(child, population);
        while (duplicate) {
            if (allowDup) {
                child = orderOneCrossover(population[winners[0]], population[winners[1]]);
            } else {
                child = orderOneCrossoverWoRepetition(population[winners[0]], population[winners[1]]);
            }
            duplicate = duplicate(child, population);

            // bunu diger case icin kullanmak istersen && kaldir
            if (duplicate && !allowDup) {
                winners = getWinners();
            }
        }

        // replace the loser with the child
        population[loser] = child;
    }

    private boolean parentLengthsOk(Integer[] p1, Integer[] p2) {
        int l = p1.length;
        if (l != p2.length) {
            System.err.println("parents must have equal lengths");
            return false;
        }

        return true;
    }

    public Integer[] orderOneCrossoverWoRepetition(Integer[] parent1, Integer[] parent2) {
        if(!parentLengthsOk(parent1, parent2)) return null;

        int l = parent1.length;
        ParentIndices parentIndices = new ParentIndices(rand.nextInt(l), rand.nextInt(l));

        // create the child .. initial elements are -1
        Integer[] child = new Integer[l];
        Arrays.fill(child, -1);
        child = Arrays.copyOfRange(parent1, parentIndices.getIndex1(), parentIndices.getIndex2());

        // array to hold elements of parent1 which are not in child yet
        int[] y = new int[l - (parentIndices.getIndex2() - parentIndices.getIndex1()) - 1];
        int j = 0;
        for (int i = 0; i < l; i++) {
            if (!searchHelp(child, parent1[i])) {
                y[j] = parent1[i];
                j++;
            }
        }

        // rotate parent2
        // number of places is the same as the number of elements after r2
        Integer[] copy = parent2.clone();
        rotate(copy, l - r2 - 1);

        // now order the elements in y according to their order in parent2
        int[] y1 = new int[l - (r2 - r1) - 1];
        j = 0;
        for (int i = 0; i < l; i++) {
            if (searchHelp(y, copy[i])) {
                y1[j] = copy[i];
                j++;
            }
        }

        // now copy the remaining elements (i.e. remaining in parent1) into
        // child
        // according to their order in parent2 .. starting after r2!
        j = 0;
        for (int i = 0; i < y1.length; i++) {
            int ci = (r2 + i + 1) % l;// current index
            child[ci] = y1[i];
        }
        return child;
    }

    private int[] orderOneCrossover(int[] parent1, int[] parent2) {
        int l = parent1.length;
        if (l != parent2.length) {
            System.err.println("parents must have equal lengths");
            return null;
        }

        // get 1 random int between 0 and size of array
        int r1 = rand.nextInt(l);

        // create the child .. initial elements are -1
        int[] child = new int[l];
        for (int i = 0; i < r1; i++) {
            child[i] = parent1[i];
        }

        for (int i = r1; i < l; i++) {
            child[i] = parent2[i];
        }
        return child;
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

    private int[] getWinners(Double winnerSimilarityRatio, int[][] population, int convergenceFactor) {
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
                winners[i] = selectParentUsingBinaryTournament(lastWinner, getSize(), chromosomePopulation.getFitnessValues());
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

    public int selectParentUsingBinaryTournament(int lastWinner, Integer popSize, Double[] popmember) {
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

    public int randomNumber(int min, int max) {
        double d = min + rand.nextDouble() * (max - min);
        return (int) d;
    }

    private class ParentIndices {
        private Integer index1;
        private Integer index2;

        ParentIndices(Integer index1, Integer index2) {
            // make sure index1 < index2
            if(index1 > index2) {
                int temp = index1;
                index1 = index2;
                index2 = temp;
            } else if (index1 > 0){
                index1--;
            } else {
                index2++;
            }

            this.index2 = index2;
            this.index1 = index1;
        }

        Integer getIndex1() { return index1; }

        Integer getIndex2() { return index2; `}
    }

    private class ChromosomeFitnessValues {
        private Double[] fitnessValues;
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

        public Double[] getFitnessValues() {
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
