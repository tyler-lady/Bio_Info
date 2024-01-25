import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class Main {

    /*
    EXPECTED OUTPUT FROM DATASET>TXT
    TCTCGGGG
    CCAAGGTG
    TACAGGCG
    TTCAGGTG
    TCCACGTG
     */

    public static void main(String[] args) throws IOException {

        String dataPath = "D:\\dontm\\Documents\\School\\B363\\hw7\\src\\dataset.txt";

        FileReader fr = new FileReader(dataPath);
        BufferedReader br = new BufferedReader(fr);

        ArrayList<String> list = new ArrayList<>();

        // Declaring a string variable
        String st;
        // Condition holds true till
        // there is character in a string
        // K, T, N, then DNA
        while ((st = br.readLine()) != null) {

            // Print the string
            list.add(st);
            System.out.println(st);
        }

        int k;
        int t;
        int n;

        String dna;

        String[] splitted = list.get(0).split("\\s+");
        System.out.println(Arrays.toString(splitted));
        k = Integer.parseInt(splitted[0]);
        t = Integer.parseInt(splitted[1]);
        n = Integer.parseInt(splitted[2]);

        ArrayList<String> DnaList = new ArrayList<>();

        for (int i = 1; i < list.size(); i++){
            DnaList.add(list.get(i));
        }
        //dna = builder.toString();

        RandomizedMotifSearch RS = new RandomizedMotifSearch(k,t,n,DnaList);
        //RS.randomizedMotifSearch();
        GibbsAligner GA = new GibbsAligner(k,t,n,DnaList);
        GA.randomizedMotifSearch();
    }
}