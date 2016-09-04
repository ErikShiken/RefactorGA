/**
 * Created by Erik on 9/4/2016.
 */
public class CrossoverParentIndices {
    private Integer index1;
    private Integer index2;

    CrossoverParentIndices(Integer index1, Integer index2) {
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

    Integer getIndex2() { return index2; }
}
