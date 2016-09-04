import java.io.File;
import java.util.*;

class GeneticAlgorithmConfiguration {
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

    GeneticAlgorithmConfiguration() {
        bestChromosomesMap = new HashMap<>();
        iterationChromosomeMap = new HashMap<>();
        allowDuplicateGenes = false;
        batchProcessing = false;
        rand = new Random();
    }

    File getFitnessFile() {
        return new File(workingPath + "\\" + FITNESS_FILE_NAME);
    }

    boolean duplicateGenesEnabled() {
        return allowDuplicateGenes;
    }

    boolean batchModeEnabled() {
        return batchProcessing;
    }

	public Iterable<String> getBestChromosomes() {
		return bestChromosomesMap.keySet();
	}
	
	Double getBestChromosomeTime(String chromosome) {
		return bestChromosomesMap.get(chromosome);
	}
	
	Map<String, Double> getBest() {
		return bestChromosomesMap;
	}
	
	Map<String, Integer> getChromosomeIterations() {
		return iterationChromosomeMap;
	}
	
	public Integer getChromosomeIteration(String chromosome) {
		return iterationChromosomeMap.get(chromosome);
	}
	
	/**
	 * @return the iterations
	 */
	int getIterations() {
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
	double getCrossover() {
		return crossover;
	}

	/**
	 * @param crossover the crossover to set
	 */
	void setCrossover(double crossover) {
		this.crossover = crossover;
	}

	/**
	 * @return the convergence
	 */
	int getConvergence() {
		return convergence;
	}

	/**
	 * @param convergence the convergence to set
	 */
	void setConvergence(int convergence) {
		this.convergence = convergence;
	}

	/**
	 * @return the similarity
	 */
	double getSimilarity() {
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
	int getSize() {
		return size;
	}

	/**
	 * @param size the size to set
	 */
	void setSize(int size) {
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
	void setInitMutation(double initMutation) {
		this.initMutation = initMutation;
		this.mutation = initMutation;
	}

	String getBestOverallChromosome() {
        if (bestChromosomeOverall == null)
            return "";

		return bestChromosomeOverall.getChromosome();
	}

	void addBestChromosome(String chromosome, Double bestTime) {
		// do not overwrite the original time seen for this chromosome
		if(bestChromosomesMap.get(chromosome) == null) {
			this.bestChromosomesMap.put(chromosome, bestTime);

            updateBestOverall(chromosome, bestTime);
        }
	}

    void updateBestOverall(String chromosome, Double bestTime) {
        if(bestChromosomeOverall == null || bestChromosomeOverall.getTime() > bestTime)
			bestChromosomeOverall = new ChromosomeTimePair(chromosome, bestTime);
	}

    void increaseMutationRate(Double v) {
        this.mutation += v;
    }

    void allowDuplicates() {
        allowDuplicateGenes = true;
    }

    void batchMode() {
        batchProcessing = true;
    }

    void setWorkingPath() {
        Scanner scan = new Scanner(System.in);
        System.out.println("Provide the path where the simulation is expecting input and provides output:");
        workingPath = scan.nextLine();
        scan.close();
    }

    String getWorkingPath() {
        return workingPath;
    }

    void setPopulationFitnessValues() {
        chromosomePopulation.updateFitnessValues(getFitnessFile());
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

    private void reproduction(Integer[] winners, int loser, boolean allowDup, Integer[][] population) {
        Integer[] child;
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
        CrossoverParentIndices parentIndices = new CrossoverParentIndices(rand.nextInt(l), rand.nextInt(l));

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

    private Integer[] orderOneCrossover(Integer[] parent1, Integer[] parent2) {
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

    private Integer[] getWinners(Double winnerSimilarityRatio, Integer[][] population, int convergenceFactor) {
        int numTournament = 2;
        // winning gene positions, the start of binary selection
        Integer[] winners = new int[numTournament];

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
            Integer[] winner1 = population[winners[0]];
            Integer[] winner2 = population[winners[1]];
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
}
