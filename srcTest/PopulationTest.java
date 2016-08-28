import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class PopulationTest {

	private static Population current;
	
	@BeforeClass
	public static void setUp() {
		current = new Population();
	}
	
	@Test
	public void testPopulation() {
		Assert.assertNotNull(current.getBest());
		Assert.assertNotNull(current.getChromosomeIterations());
	}

	@Test
	public void testAddBestChromosome() {
		String chromosome = "1,2,3,4,5,6,7,8,9";
		Double bestTime = 200d;
		current.addBestChromosome(chromosome, bestTime);
		current.addBestChromosome(chromosome, bestTime * 2);
		
		Assert.assertEquals(current.getBestChromosomeTime(chromosome), bestTime);
	}
	
	@Test
	public void testGetBestOverallChromosome() {
		String chromosome = "1,2,3,4,5,6,7,8,9";
		Double bestTime = 200d;
		current.addBestChromosome(chromosome, bestTime);
		
		Assert.assertTrue(current.getBestOverallChromosome().equalsIgnoreCase(chromosome));
	}

    @Test
    public void testUpdateBestOverallCreateNewBest() {
        String chromosome = "1,2,3,4,5,6,7,8,9";
        Double bestTime = 200d;

        current.updateBestOverall(chromosome, bestTime);

        Assert.assertEquals("Current best mismatch", current.getBestOverallChromosome(), chromosome);
    }
}