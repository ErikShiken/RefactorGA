
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;

public class RefactoredGA {
	private static Population currentPopulation;
	// fitness value file
	private static int StoppingCondition = 0;
	
	private static File populationOutput;
	private static Scanner scan;
	private static File fitFile;

	public static int fitnessValue = 0;

	static Random rand = new Random();

	public static void main(String[] args) throws Exception {
		currentPopulation = new Population();
		// bestChromosomesMap = new HashMap<String, Double>();
		// iterationChromosomeMap = new HashMap<>();
		// boolean allowDupWorkers = true;
		// runGA(allowDupWorkers);
		// printResults();

		File input = new File("settings.txt");
		Scanner inputScanner = new Scanner(input);
		String currentLine;
		
		Scanner lineParser = null;
		while(inputScanner.hasNext()){
			currentLine = inputScanner.nextLine();
			lineParser = new Scanner(currentLine);
			
			currentPopulation.setSize(Integer.parseInt(lineParser.next()));
			currentPopulation.setCrossover(Float.parseFloat(lineParser.next()));
			currentPopulation.setInitMutation(Float.parseFloat(lineParser.next()));
			
			boolean allowDupWorkers = true;
			boolean batch = false;
			runGA(allowDupWorkers, batch);
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
		for (String chromosome : currentPopulation.getBestChromosomes()) {
			results.append(String.format("%s : ", chromosome));
			results.append(String.format("%f : ", currentPopulation.getBestChromosomeTime(chromosome)));
			results.append(String.format("%d", currentPopulation.getChromosomeIteration(chromosome)));
			results.append(System.getProperty("line.separator"));
		}
		
		results.append(String.format("XORate = %f", currentPopulation.getCrossover()));
		results.append(System.getProperty("line.separator"));
		
		results.append(String.format("mutRate = %f", currentPopulation.getMutation()));
		results.append(System.getProperty("line.separator"));
		
		return results.toString();
	}
	
	private static void printToFile(File output) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		
		String bestChromosomeResults = getBestResults();

		String best = currentPopulation.getBestOverallChromosome();
				
		bw.write(bestChromosomeResults);
		bw.write(String.format("Best chromosome : %s : Iteration = %d : Score = %f : mutRate = %f", 
				best, currentPopulation.getChromosomeIteration(best), 
				currentPopulation.getBestChromosomeTime(best), 
				currentPopulation.getMutation()));
		bw.newLine();
		bw.flush();
		bw.close();

	}
	
	private static void printResultsToFile() throws IOException {
		File newOutput = 
			createResultsFile(
				currentPopulation.getSize(), 
				currentPopulation.getCrossover(),
				currentPopulation.getInitMutation());
		
		printToFile(newOutput);
	}

	public static void runGA(boolean allowDupGenes, boolean batch) throws Exception {
		if(allowDupGenes) System.out.println("Allowing duplicates");
		else System.out.println("Duplicates NOT allowed.");
		String path = "C:\\Users\\etest\\Desktop\\AirInterdiction\\";
		runSetup(path, allowDupGenes, batch);

		// convergentRun keeps track of how many iterations have occurred where a new
		// solution has not been found
		int convergentRun = 0;
		for (int loop = 0; loop < currentPopulation.getIterations(); loop++) {
			if(loop % 100 == 0) System.out.println("loop = " + loop);
			
			double largestFitness = Double.MIN_VALUE;
			double smallestFitness = Double.MAX_VALUE;
			int loser = -1;
			int winner = -1;

			Scanner fitnessScanner = new Scanner(fitFile);
            for (int i = 0; i < currentPopulation.getSize(); i++) {
                popmember[i] = fitnessScanner.nextDouble();
				if (popmember[i] > largestFitness) {
					largestFitness = popmember[i];
					loser = i;
				} 
				
				if (popmember[i] < smallestFitness) {
					smallestFitness = popmember[i];
					winner = i;
				}
			}
			fitnessScanner.close();

			if(updateWinner(winner, loop))
			{	
				System.out.println("convergentRun was " + convergentRun);
				convergentRun = 0;
			}
			else {	convergentRun++; }
			
			// stopping condition
			if (popmember[winner] < StoppingCondition) {
				System.out.println("Stopping condition reached");
				break;
            } else if (convergentRun >= currentPopulation.getConvergence()) {
                System.out.println("GA has converged based on the given convergenceFactor ("
                        + currentPopulation.getConvergence() + ")");
                break;
            } else if (convergentRun == currentPopulation.getConvergence() / 2) {
                currentPopulation.increaseMutationRate(0.01);
                System.out.println("mutRate updated : " + currentPopulation.getConvergence());
            }

            for (int i = 0; i < currentPopulation.getSize(); i++) {
                if (rand.nextDouble() < currentPopulation.getCrossover()) {
                    reproduction(getWinners(), loser, allowDupGenes);
				}
			}

			// Mutate the new population a bit to add some new genetic material
			mutate(winner);

			if (batch) {
                writeNewPopulationToFile();

				runSimulation(path);
			} else {
				writeNewPopulationIterative(path);
			}
		}

        if (convergentRun < currentPopulation.getConvergence())
            System.out.println("GA simulation completed after " + currentPopulation.getIterations() + " rounds.");
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

    private static int[] getWinners(Double winnerSimilarityRatio, int[][] population, int convergenceFactor) {
        int numTournament = 2;
		// winning gene positions, the start of binary selection
		int[] winners = new int[numTournament];

        int tempConverge = convergenceFactor;
        double diffCount = 0;
		double tempSimilarity = winnerSimilarityRatio;
		while (diffCount / population[0].length <= tempSimilarity) {
			if(diffCount > 0) {
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

    private static void reproduction(int[] winners, int loser, boolean allowDup, int[][] population) {
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

	// if this winner has not been seen, add it to the bestChromosomeMap
	// -- returns true if the map has been updated, false otherwise
    private static boolean updateWinner(int winner, int loop, int[][] populaton) {
        boolean updated = false;
		String chromString = convertToString(population[winner]);
		if (bestChromosomesMap.get(chromString) == null) {
			System.out.println("adding to map : " + chromString + "(" + popmember[winner] + ") : loop = " + loop);
			bestChromosomesMap.put(chromString, popmember[winner]);
			iterationChromosomeMap.put(chromString, (loop + 1));
			updated = true;
		}
		
		return updated;
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

	private static boolean duplicate(int[] child, int[][] population) {
		boolean foundDup = true;
		for (int i = 0; i < population.length; i++) {
			foundDup = true;
			for (int j = 0; j < child.length; j++) {
				if (population[i][j] != child[j]) {
					foundDup = false;
					break;
				}
			}

			// if there is a duplicate, break out of the loop and return true.
			// no more work is needed
			if (foundDup)
				break;
		}
		return foundDup;
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

	private static int[] orderOneCrossover(int[] parent1, int[] parent2) {
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

	public static int[] orderOneCrossoverWoRepetition(int[] parent1, int[] parent2) {
		int l = parent1.length;
		if (l != parent2.length) {
			System.err.println("parents must have equal lengths");
			return null;
		}
		// get 2 random ints between 0 and size of array
		int r1 = rand.nextInt(l);
		int r2 = rand.nextInt(l);

		// to make sure the r1 < r2
		while (r1 >= r2) {
			r1 = rand.nextInt(l);
			r2 = rand.nextInt(l);
		}

		// create the child .. initial elements are -1
		int[] child = new int[l];
		for (int i = 0; i < l; i++) {
			child[i] = -1;
		}

		// copy elements between r1, r2 from parent1 into child
		for (int i = r1; i <= r2; i++) {
			child[i] = parent1[i];
		}

		// array to hold elements of parent1 which are not in child yet
		int[] y = new int[l - (r2 - r1) - 1];
		int j = 0;
		for (int i = 0; i < l; i++) {
			if (!searchHelp(child, parent1[i])) {
				y[j] = parent1[i];
				j++;
			}
		}

		// rotate parent2
		// number of places is the same as the number of elements after r2
		int[] copy = parent2.clone();
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

	public static void rotate(int[] arr, int order) {
		int offset = arr.length - order % arr.length;
		if (offset > 0) {
			int[] copy = arr.clone();
			for (int i = 0; i < arr.length; ++i) {
				int j = (i + offset) % arr.length;
				arr[i] = copy[j];
			}
		}
	}

	private static boolean searchHelp(int[] arr, int key) {
		for (int i : arr) {
			if (i == key)
				return true;
		}

		return false;
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
}// End RefactorGA class