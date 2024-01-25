import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static java.lang.Long.sum;

public class GibbsAligner {
    int k,t,n;
    ArrayList<String> dna;

    GibbsAligner(int k, int t, int n, ArrayList<String> dna){
        this.k = k;
        this.t = t;
        this.n = n;
        this.dna = dna;
    }

    /*
    GIBBSSAMPLER(Dna, k, t, N)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs
                       except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
     */

    /*
    EXPECTED OUTPUT FROM DATASET>TXT
    TCTCGGGG
    CCAAGGTG
    TACAGGCG
    TTCAGGTG
    TCCACGTG
    //917606
     */

    public void randomizedMotifSearch() throws IOException {
        //2nd, randomly select T_seq K_mers from DNA_List;
        ArrayList<String> motifs = randomSelectMotifs(dna, k, t);
        //ArrayList<String> bestMotifs = new ArrayList<String>(motifs);

        //get the motifsScore of bestMotifs;
        //int bestScore = motifsScore(motifs);
        ArrayList<String> current_motifs = new ArrayList<>();
        ArrayList<String> bestMotifs = new ArrayList<String>(current_motifs);

        for(int i=0;i<20;i++){
            current_motifs = updateMotifs(motifs, dna, k, t);
            if(i == 0 || motifsScore(current_motifs) < motifsScore(bestMotifs))
                bestMotifs = current_motifs;
        }

        //3rd, update motifs and bestMotifs
        //ArrayList<String> bestMotifs = updateMotifs(motifs, dna, k, t);


        //4th, printout bestMotifs arrayList:
        printArrayList(bestMotifs);


        String pathname = "D:\\dontm\\Documents\\School\\B363\\hw7\\src\\";
        //5th, output writer into an randomMotifSearch.txt document, with carriage return
        String doc_out = "randomMotifSearch.txt";
        File output_file = new File(pathname + doc_out);
        BufferedWriter output = new BufferedWriter(new FileWriter(output_file));

        for(int i=0; i<bestMotifs.size(); i++){
            output.write(bestMotifs.get(i) + "\r");
        }

        //close writer
        output.close();
        //end main();
    }


    /*********************************
     * Printout an arrayList of strings;
     * @param motifs
     */
    private static void printArrayList(ArrayList<String> motifs) {
        System.out.println("\n Printout the Motifs arrayList:");

        for(int i=0; i<motifs.size(); i++){

            System.out.print(motifs.get(i) + "\r\n");

        }//end for i<motifs.size() loop;

    }//end of pintArrayList() method;

    private String profileRandomizedKmer(String dna, int k, double[][] profile){
        /*
        nuc_loc = {nucleotide: index for index, nucleotide in enumerate('ACGT')}
    probs = []
    for i in range(len(dna) - k):
        current_prob = 1.
        for j, nucleotide in enumerate(dna[i:i + k]):
            current_prob *= prof[j][nuc_loc[nucleotide]]
        probs.append(current_prob)

    i = numpy.random.choice(len(probs), p = numpy.array(probs) / numpy.sum(probs))
    return dna[i:i + k]
         */
        HashMap<Character, Integer> nuc = new HashMap<>();
        nuc.put('A', 0);
        nuc.put('C', 1);
        nuc.put('G', 2);
        nuc.put('T', 3);

        ArrayList<Double> probs = new ArrayList<>();
        for(int i=0; i<dna.length()-k;i++){
            double currProb = 1;

            for(int j=0;j<4;j++){
                for(Character nucleotide: dna.substring(i,i+k).toCharArray()){
                    currProb *= profile[j][nuc.get(nucleotide)];
                }

            }
            probs.add(currProb);
        }
        RandomCollection rc = new RandomCollection<>();
        for(int i = 0; i < probs.size(); i++){
            rc.add((probs.get(i)/sumList(probs))*100, i);
        }

        int index = (int) rc.next();

        return dna.substring(index, index+k);
    }

    private double sumList(ArrayList<Double> list){
        double total = 0;
        for(double d : list){
            total += d;
        }
        return total;
    }

    private ArrayList<String> updateMotifs(ArrayList<String> motifs, ArrayList<String> DNAlist, int k_mer, int t_seq) {
        motifs = randomSelectMotifs(DNAlist, k_mer, t_seq);
        ArrayList<String> bestMotifs = new ArrayList<String>(motifs);
        int bestScore = motifsScore(motifs);

        for(int k = 1; k<this.n;k++) {
            //motifs = randomSelectMotifs(DNAlist, k_mer, t_seq);
            int i = ThreadLocalRandom.current().nextInt(0, t_seq - 1);

            double[][] profile = createProMatrix(motifs);
            //get a new motifs matrix according to the profileMatrix we just got
            ArrayList<String> currMotifs = new ArrayList<String>(motifs);
            System.out.println("i is: " + i);
            currMotifs.set(i, profileRandomizedKmer(DNAlist.get(i), k_mer, profile));

//                for(int j=0; j<t_seq; j++){
//                    int i = ThreadLocalRandom.current().nextInt(0, t_seq-1);
//                    System.out.println("i is: " + i);
//                    currMotifs.set(i, profileRandomizedKmer(DNAlist.get(i), k_mer, profile));
//                }
        }
                printArrayList(currMotifs);
                int currScore = motifsScore(currMotifs);
                if(currScore > bestScore){
                    bestScore = currScore;
                    motifs = currMotifs;
                    bestMotifs = new ArrayList<String>(motifs);
                }


        return bestMotifs;
    }

    /******************************
     * motif matrix formed by first k-mers in each string from DNA sequences;
     * @param dNA_List
     * @param k_mer
     * @return
     */
    private static ArrayList<String> generateBestMotifs(ArrayList<String> DNA_list, int k_mer) {
        // TODO Auto-generated method stub
        ArrayList<String> kmerList = new ArrayList<String>();

        for(int i=0; i<DNA_list.size(); i++){

            String subSeq = DNA_list.get(i).substring(0, k_mer);
            kmerList.add(subSeq);

        }//end for i<DNA_list.size() loop;

        return kmerList;

    }//end generateBestMotifs() method;

    /******************************
     * get a best match sub-sequence from the profileMatrix and the given original sequence;
     * @param profileMatrix
     * @param sequence
     * @param k_mer
     * @return
     */
    private static String profileMostProbable(double[][] profileMatrix,	String sequence, int k_mer) {
        /************
         * profileMatrix
         * A
         * T
         * G
         * C
         */

        //1st, get 1st k_mer from this sequence
        String Kmer = sequence.substring(0, k_mer);
        double score = profileScore(Kmer, profileMatrix);

        int Len = sequence.length();
        int lastStartPoint = Len - k_mer + 1;

        for(int i=1; i<lastStartPoint; i++){

            String currSeq = sequence.substring(i, i+k_mer);

            double currPro = profileScore(currSeq, profileMatrix);
            if(currPro > score){

                Kmer = currSeq;
                score = currPro;
            }

        }//end for i<lastStartPoint loop;

        return Kmer;

    }

    private static double profileScore(String seq, double[][] profileMatrix) {
        /************
         * profileMatrix
         * A
         * T
         * G
         * C
         */
        double score = 1;
        int Len = seq.length();

        for(int i=0; i<Len; i++){

            double currPro = 0;
            switch(seq.charAt(i)){

                case 'A': currPro = profileMatrix[0][i]; break;
                case 'T': currPro = profileMatrix[1][i]; break;
                case 'G': currPro = profileMatrix[2][i]; break;
                case 'C': currPro = profileMatrix[3][i]; break;

            }//end switch-case loop;

            score = score * currPro;

        }//end for i<Len loop;

        return score;

    }//end profileScore() method;

    private static double[][] createProMatrix(ArrayList<String> motifs) {
        // TODO Auto-generated method stub
        int row = motifs.size();
        int col = motifs.get(0).length();

        /*****
         * A =-=-=-=
         * T =-=-=-=
         * G =-=-=-=
         * C =-=-=-=
         */
        double[][] profile = new double[4][col];

        for (int i = 0; i < col; i++) {
            double total = motifs.size();
            double A = 1;
            double T = 1;
            double G = 1;
            double C = 1;

            for (int j = 0; j < row; j++) {

                switch (motifs.get(j).charAt(i)) {

                    case 'A':
                        A++;
                        break;
                    case 'T':
                        T++;
                        break;
                    case 'G':
                        G++;
                        break;
                    case 'C':
                        C++;
                        break;

                }

            }//end j<row loop;

            profile[0][i] = A / total;
            profile[1][i] = T / total;
            profile[2][i] = G / total;
            profile[3][i] = C / total;

        }//end for i<col loop;

        return profile;
    }

    /*********************************
     * get a score of one motif matrix
     * @param motifList
     * @return
     */
    private static int motifsScore(ArrayList<String> motifList) {
        int row = motifList.size();
        int col = motifList.get(0).length();
        int score = 0;
        System.out.println("Cols" + col);
        for(int i=0; i<col; i++){

            int A = 0;
            int T = 0;
            int G = 0;
            int C = 0;

            for(int j=0; j<row; j++){
                switch(motifList.get(j).charAt(i)){
                    case 'A': A++; break;
                    case 'T': T++; break;
                    case 'G': G++; break;
                    case 'C': C++; break;
                }//end switch-case;

            }//end for j<row loop;


            int max = A;
            if( T > max) max = T;
            if( G > max) max = G;
            if( C > max) max = C;

            //System.out.println("A: " + A +", T: " + T + ", G: " + G +", C: " + C  +", max: " + max);

            score = score + row - max;

        }//end for i<col loop;

        return score;

    }//end of motifsScore() method;

    /*********************
     * randomly select a k_mer sub-sequence from each DNA sequence;
     * put them into an arrayList: motifs
     * @param DNAlist
     * @param k_mer
     * @param t_seq
     * @return motifs
     */
    private static ArrayList<String> randomSelectMotifs(ArrayList<String> DNAlist, int k_mer, int t_seq) {
        ArrayList<String> motifs = new ArrayList<String>();

        int Len = DNAlist.get(0).length();

        int start = 0;
        int end = Len - k_mer;
        //	System.out.println("len: " + Len + ", start: " + start +", end: " + end);

        for(int i=0; i<t_seq; i++){

            int rand = (int) (Math.random() * (end-start));

            String subSeq = DNAlist.get(i).substring(rand, rand + k_mer);
            motifs.add(subSeq);

            //	System.out.println("Rand:\t" + rand + "\t  Seq: " + subSeq);
        }

        return motifs;

    }//end randomSelectMotifs() method;

}
