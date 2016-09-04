import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Scanner;

public class RefactoredGA {
    // fitness value file
    private static Double StoppingCondition = 0d;
	private static Scanner scan;

	public static void main(String[] args) throws Exception {
        GeneticAlgorithmConfiguration config = new GeneticAlgorithmConfiguration();
        config.allowDuplicates();
        config.batchMode();

		File input = new File("settings.txt");
		Scanner inputScanner = new Scanner(input);
		String currentLine;
		
		Scanner lineParser = null;
		while(inputScanner.hasNext()){
			currentLine = inputScanner.nextLine();
			lineParser = new Scanner(currentLine);

            config.setSize(Integer.parseInt(lineParser.next()));
            config.setCrossover(Float.parseFloat(lineParser.next()));
            config.setInitMutation(Float.parseFloat(lineParser.next()));
            config.setWorkingPath();

            runGA(config);
            printResultsToFile();
			
			if(scan != null)
				scan.close();
		}
		
		if(lineParser != null)
			lineParser.close();
		if(inputScanner != null)
			inputScanner.close();
	}

	private static File createResultsFile(int popSize, double XORate, double initMutRate) {
		Calendar cal = Calendar.getInstance();
		StringBuilder output = new StringBuilder(popSize);
		output.append(String.format("_%.2f_", XORate));
		output.append(String.format("%.2f_", initMutRate));
		
		String dfn = "yy.MM.dd.HH.mm.ss";
		SimpleDateFormat sdf = new SimpleDateFormat(dfn);
		output.append(sdf.format(cal.getTime()));
		
		return new File(output.toString());
	}
	
	private static String getBestResults() {
		StringBuffer results = new StringBuffer();
        for (String chromosome : config.getBestChromosomes()) {
            results.append(String.format("%s : ", chromosome));
            results.append(String.format("%f : ", config.getBestChromosomeTime(chromosome)));
            results.append(String.format("%d", config.getChromosomeIteration(chromosome)));
            results.append(System.getProperty("line.separator"));
        }

        results.append(String.format("XORate = %f", config.getCrossover()));
        results.append(System.getProperty("line.separator"));

        results.append(String.format("mutRate = %f", config.getMutation()));
        results.append(System.getProperty("line.separator"));

        return results.toString();
	}
	
	private static void printToFile(File output) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		
		String bestChromosomeResults = getBestResults();

        String best = config.getBestOverallChromosome();

        bw.write(bestChromosomeResults);
        bw.write(String.format("Best chromosome : %s : Iteration = %d : Score = %f : mutRate = %f",
                best, config.getChromosomeIteration(best),
                config.getBestChromosomeTime(best),
                config.getMutation()));
        bw.newLine();
        bw.flush();
		bw.close();

	}
	
	private static void printResultsToFile() throws IOException {
		File newOutput = 
			createResultsFile(
                    config.getSize(),
                    config.getCrossover(),
                    config.getInitMutation());

        printToFile(newOutput);
    }

    public static void runGA(GeneticAlgorithmConfiguration config) throws Exception {
        if (config.duplicateGenesEnabled()) System.out.println("Allowing duplicates");
        else System.out.println("Duplicates NOT allowed.");

        runSetup(config.getWorkingPath(), config.duplicateGenesEnabled(), config.batchModeEnabled(), config.getSize());

		// convergentRun keeps track of how many iterations have occurred where a new
		// solution has not been found
		int convergentRun = 0;
        for (int loop = 0; loop < config.getIterations(); loop++) {
            if (loop % 100 == 0) System.out.println("loop = " + loop);

            config.setPopulationFitnessValues();

			if(config.updateWinner(loop))
			{	
				System.out.println("convergentRun was " + convergentRun);
				convergentRun = 0;
			}
			else {	convergentRun++; }
			
			// stopping condition
			if (config.getWinningFitnessValue() < StoppingCondition) {
				System.out.println("Stopping condition reached");
				break;
            } else if (convergentRun >= config.getConvergence()) {
                System.out.println("GA has converged based on the given convergenceFactor ("
                        + config.getConvergence() + ")");
                break;
            } else if (convergentRun == config.getConvergence() / 2) {
                config.increaseMutationRate(0.01);
                System.out.println("mutRate updated : " + config.getConvergence());
            }

            config.doCrossover();

			config.mutate();
			// Mutate the new population a bit to add some new genetic material
			mutate(winner);

			if (config.batchModeEnabled()) {
                writeNewPopulationToFile(config.getWorkingPath());

				runSimulation(config.getWorkingPath());
			} else {
				writeNewPopulationIterative(config.getWorkingPath(), loop);
			}
		}

        if (convergentRun < config.getConvergence())
            System.out.println("GA simulation completed after " + config.getIterations() + " rounds.");
        else
			System.out.println("converged");
	}

    private static void writeNewPopulationIterative(String path, Integer popSize) throws Exception {
        float[] averages = new float[popSize];

		// for entire population
		for (int i = 0; i < popSize; i++) {
			writeGeneToFile(i);

			runSimulation(path);
			
			updateAverages(averages, i);
			
			if(i % 2 == 0) System.out.print(".");
		}
		
		System.out.println();
		fitFile.delete();
		fitFile.createNewFile();
		FileWriter averagesWriter = new FileWriter(fitFile);
		BufferedWriter bufferedAverages = new BufferedWriter(averagesWriter);
		for (float f : averages) {
			bufferedAverages.write(Float.toString(f));
			bufferedAverages.newLine();
		}
		bufferedAverages.close();
	}

	private static void updateAverages(float[] averages, int index) throws FileNotFoundException {
		scan = new Scanner(fitFile);
		float sum = 0f;
		int numLines = 0;
		while (scan.hasNextLine()) {
			sum += Float.parseFloat(scan.nextLine());
			numLines++;
		}

		averages[index] = sum / numLines;
		scan.close();
	}

    private static void runSetup(String path, boolean allowDupGenes, boolean batch, Integer popSize) throws Exception {
        if (batch) {
			// write population to a file
			if (allowDupGenes) {
				initPopulation();
			} else {
				initPopulationWoDuplicates();
			}

			// run the simulation
			runSimulation(path);
		} else {
			if (allowDupGenes) {
				initIterativePopulation();
			} else {
				initIterativePopulationWoDuplicates();
			}

			float[] averages = new float[popSize];
			for (int i = 0; i < popSize; i++) {
				writeGeneToFile(i);

				runSimulation(path);

				updateAverages(averages, i);
				
				if (i % 2 == 0) System.out.print(".");
			}

			System.out.println();
			fitFile.delete();
			fitFile.createNewFile();
			FileWriter averagesWriter = new FileWriter(fitFile);
			BufferedWriter bufferedAverages = new BufferedWriter(averagesWriter);
			for (float f : averages) {
				bufferedAverages.write(Float.toString(f));
				bufferedAverages.newLine();
			}
			bufferedAverages.close();
		}
	}

    private static void initIterativePopulationWoDuplicates(Integer numAgents, Integer popSize, int[][] population) {
        // Agents (must be equal to numAgents)
		int[] Agents = new int[numAgents];
		for (int i = 0; i < numAgents; i++) {
			Agents[i] = i + 1;
		}

		// for entire population
		for (int i = 0; i < popSize; i++) {
			shuffle(Agents);

			boolean foundDuplicate = duplicate(Agents, population);
			while (foundDuplicate) {
				Agents = random(Agents.length);
				foundDuplicate = duplicate(Agents, population);
			}

			for (int j = 0; j < numAgents; j++) {
				population[i][j] = Agents[j];
			}
		}
	}

    private static void writeGeneToFile(int index, int[][] population) throws Exception {
        populationOutput.delete();
		populationOutput.createNewFile();
		StringBuilder sb = new StringBuilder("");
		FileWriter fw = new FileWriter(populationOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		
		for (int i = 0; i < population[index].length; i++) {
			sb.append(population[index][i] + " ");
		}
		
		bw.write(sb.toString());
		bw.write(System.getProperty("line.separator"));
		sb.setLength(0);
		bw.flush();
		bw.close();
	}

    private static void initIterativePopulation(Integer numAgents, Integer popSize, int[][] population) {
        // Agents (must be equal to numAgents)
		int[] Agents = new int[numAgents];
		for (int i = 0; i < numAgents; i++) {
			Agents[i] = i + 1;
		}

		// for entire population
		for (int i = 0; i < popSize; i++) {
			shuffle(Agents);

			for (int j = 0; j < numAgents; j++) {
				population[i][j] = Agents[j];
			}
		}
	}

    private static void mutate(int winner, Integer popSize, int[][] population, Double mutRate) {
        int[] tempChromosome = null;
		StringBuilder winnerStr = new StringBuilder();
		for (int i = 0; i < popSize; i++) {
			// skip this population member if it is one of the winners
			if(winner == i) {
				for(int j = 0; j < population[i].length; j++ ) {
					winnerStr.append(population[i][j]);
					if( j < population[i].length - 1)
						winnerStr.append(",");
				}
				continue;
			}
			
			if (rand.nextDouble() < mutRate) {
				boolean duplicate = true;
				while (duplicate) {
					tempChromosome = swapMutation(population[i]);
					if (!duplicate(tempChromosome, population))
						duplicate = false;

				}
				population[i] = tempChromosome;
			}
		}
	}

	private static void runSimulation(String path) throws IOException, InterruptedException {

	}

	private static String convertToString(int[] chromosome) {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < chromosome.length; i++) {
			sb.append(chromosome[i]);
			if (i + 1 < chromosome.length)
				sb.append(",");
		}

		return sb.toString();
	}

    private static void writeNewPopulationToFile(int[][] population, Integer popSize, Integer numAgents)
            throws IOException {
        populationOutput.delete();
		populationOutput.createNewFile();
		StringBuilder sb = new StringBuilder("");
		FileWriter fw = new FileWriter(populationOutput);
		BufferedWriter bw = new BufferedWriter(fw);

		// for entire population
		for (int i = 0; i < popSize; i++) {
			for (int j = 0; j < numAgents; j++)
				sb.append(population[i][j] + " ");

			bw.write(sb.toString());
			bw.write(System.getProperty("line.separator"));
			sb.setLength(0);
		}
		bw.close();
	}

	/**
	 * performs swap mutation on an array of ints
	 *
	 * @param parent
	 *            the int array
	 * @return array the mutated array
	 */
	private static int[] swapMutation(int[] parent) {
		int[] array = parent.clone();
		int l = array.length;
		// get 2 random integers between 0 and size of array
		int r1 = rand.nextInt(l);
		int r2 = rand.nextInt(l);
		// to make sure the 2 numbers are different
		while (r1 == r2)
			r2 = rand.nextInt(l);

		// swap array elements at those indices
		int temp = array[r1];
		array[r1] = array[r2];
		array[r2] = temp;

		return array;
	}

	// initialize the population and allow duplicate genes in a chromosome
    private static void initPopulationWoDuplicates(Integer numAgents, Integer popSize, int[][] population) throws IOException {
        // Agents (must be equal to numAgents)
		int[] Agents = new int[numAgents];
		for (int i = 0; i < numAgents; i++) {
			Agents[i] = i + 1;
		}

		StringBuilder sb = new StringBuilder("");
		FileWriter fw = new FileWriter(populationOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		// for entire population
		for (int i = 0; i < popSize; i++) {
			shuffle(Agents);

			for (int j = 0; j < numAgents; j++) {
				population[i][j] = Agents[j];
				sb.append(population[i][j] + " ");
			}

			bw.write(sb.toString());
			bw.write(System.getProperty("line.separator"));
			sb.setLength(0);
		}
		bw.close();

	}

	// create an initial gene without duplicate chromosomes
    private static void initPopulation(Integer popSize, int[][] population, Integer numAgents) throws IOException {
        // Agents (must be equal to numAgents)
		int[] Agents = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

		StringBuilder sb = new StringBuilder("");
		FileWriter fw = new FileWriter(populationOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		// for entire population
		for (int i = 0; i < popSize; i++) {
			boolean foundDuplicate = true;
			while (foundDuplicate) {
				Agents = random(Agents.length);
				foundDuplicate = duplicate(Agents, population);
			}

			population[i] = Agents;

			for (int j = 0; j < numAgents; j++) {
				sb.append(population[i][j] + " ");
			}

			bw.write(sb.toString());
			bw.write(System.getProperty("line.separator"));
			sb.setLength(0);
		}
		bw.close();

	}

	private static int[] random(int numAgents) {
		int[] array = new int[numAgents];
		for (int i = 0; i < numAgents; i++) {
			int random = rand.nextInt(numAgents) + 1;
			array[i] = random;
		}

		return array;
	}

	static void shuffle(int[] array) {
		int n = array.length;
		for (int i = 0; i < array.length; i++) {
			// Get a random index of the array past i.
			int random = i + (int) (Math.random() * (n - i));
			// Swap the random element with the present element.
			int randomElement = array[random];
			array[random] = array[i];
			array[i] = randomElement;
		}
	}
}// End RefactorGA class