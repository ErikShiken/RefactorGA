import java.util.HashMap;
import java.util.Map;

public class Population {
	private int iterations;
	private double mutation;
	private double crossover;
	private int convergence;
	private double similarity;
	private double[] popmember;
	private int[][] population;
	private int numAgents;
	private int size;
	
	// key = chromosome, value = time
	private Map<String, Double> bestChromosomesMap;
	// key = chromosome, value = iteration
	private Map<String, Integer> iterationChromosomeMap;
	
	public Population() {
		bestChromosomesMap = new HashMap<>();
		iterationChromosomeMap = new HashMap<>();
	}

	public Map<String, Double> getBest() {
		return bestChromosomesMap;
	}
	
	public Map<String, Integer> getChromosomeIterations() {
		return iterationChromosomeMap;
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
	 * @return the popmember
	 */
	public double[] getPopmember() {
		return popmember;
	}

	/**
	 * @param popmember the popmember to set
	 */
	public void setPopmember(double[] popmember) {
		this.popmember = popmember;
	}

	/**
	 * @return the population
	 */
	public int[][] getPopulation() {
		return population;
	}

	/**
	 * @param population the population to set
	 */
	public void setPopulation(int[][] population) {
		this.population = population;
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
}
