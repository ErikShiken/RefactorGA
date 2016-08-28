import java.util.HashMap;
import java.util.Map;

public class GeneticAlgorithmConfiguration {
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

    public GeneticAlgorithmConfiguration() {
        bestChromosomesMap = new HashMap<>();
        iterationChromosomeMap = new HashMap<>();
        chromosomePopulation = new ChromosomeFitnessValues();
        allowDuplicateGenes = false;
        batchProcessing = false;
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

    private class ChromosomeFitnessValues {
        private double[] populationMembers;
        private int[][] population;

        public double[] getPopulationMembers() {
            return populationMembers;
        }

        public void setPopulationMembers(double[] populationMembers) {
            this.populationMembers = populationMembers;
        }

        public int[][] getPopulation() {
            return population;
        }

        public void setPopulation(int[][] population) {
            this.population = population;
        }
    }
}
