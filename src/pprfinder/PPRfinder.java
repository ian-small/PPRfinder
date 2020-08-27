/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 *
 * @author ian
 */
public class PPRfinder {
    static boolean VERBOSE = false;

    public static HashMap<String, PPR> PPRs = new HashMap<String, PPR>();
    public static HashMap<String, Boolean> Gapless_Motif_Combinations = new HashMap<String, Boolean>();
    public static HashMap<String, Boolean> Motif_Combinations = new HashMap<String, Boolean>();

    /**
     * @param args the command line arguments args[0] protein sequence file (fasta)
     *             args[1] --domtout table from HMMER3.2
     */

    public static void main(String[] args) {

        // check arguments
        if (args.length != 2 || !(new File(args[0]).exists()) || !(new File(args[1]).exists())) {
            System.out.println("Program: PPRfinder (for characterising PPR motifs in proteins)");
            System.out.println("Version: 1.0");
            System.out.println("Usage: java -jar PPRfinder.jar <orfs.fasta> <orfs.domt>");
            System.out.println("<orfs.fasta>: fasta file of protein sequences");
            System.out.println("<orfs.domt>: domain table output from running hmmsearch on <orfs.fasta> with PPR HMMs");
            System.exit(0);
        }

        try {
            // read --domtout table
            if (PPRfinder.VERBOSE)
                System.out.println("reading " + args[1]);
            BufferedReader reader = new BufferedReader(new FileReader(args[1]));
            String s;
            PPR p;
            String[] fields;

            while ((s = reader.readLine()) != null) {
                if (s.startsWith("#")) {
                    continue;
                }
                fields = s.split(" +");
                p = PPRs.get(fields[0]);
                if (p == null) {
                    p = new PPR(fields[0]);
                    PPRs.put(fields[0], p);
                }
                p.addMotif(new Motif(fields));
            }
            reader.close();
            if (PPRfinder.VERBOSE)
                System.out.println("read " + PPRs.size() + " PPRs");

            // read fasta sequence file
            if (PPRfinder.VERBOSE)
                System.out.println("reading sequences from " + args[0]);
            reader = new BufferedReader(new FileReader(args[0]));
            StringBuilder sb = new StringBuilder(3000);
            p = null;
            while ((s = reader.readLine()) != null) {
                if (s.startsWith(">")) {
                    // finish previous PPR if there is one
                    if (p != null) {
                        p.setSequence(sb.toString());
                        sb.setLength(0);
                    }
                    fields = s.substring(1).split(" ");
                    p = PPRs.get(fields[0]);
                } else if (p != null) {
                    sb.append(s);
                }
            }
            if (p != null) {
                p.setSequence(sb.toString());
                sb.setLength(0);
            }
            reader.close();

            // read valid motif combinations
            InputStream is = PPRfinder.class.getClassLoader().getResourceAsStream("resources/motif_combinations.txt");
            reader = new BufferedReader(new InputStreamReader(is));
            while ((s = reader.readLine()) != null) {
                fields = s.split("\t");
                Motif_Combinations.put(fields[0], fields[1].equals("1") ? Boolean.TRUE : Boolean.FALSE);
            }
            reader.close();

            is = PPRfinder.class.getClassLoader().getResourceAsStream("resources/gapless_motif_combinations.txt");
            reader = new BufferedReader(new InputStreamReader(is));
            while ((s = reader.readLine()) != null) {
                fields = s.split("\t");
                Gapless_Motif_Combinations.put(fields[0], fields[1].equals("1") ? Boolean.TRUE : Boolean.FALSE);
            }
            reader.close();

            // System.out.println("inferring motifs... ");
            BufferedWriter motif_writer = new BufferedWriter(new FileWriter(args[1] + "_motifs.txt"));
            BufferedWriter orf_writer = new BufferedWriter(new FileWriter(args[1] + "_pprs.fa"));
            BufferedWriter beads_writer = new BufferedWriter(new FileWriter(args[1] + "_beads.txt"));
            // Blosum62 blosum = new Blosum62();
            // blosum.random_associations();

            // String beads;
            double score;

            for (PPR ppr : PPRs.values()) {
                if (ppr.sequence != null && !ppr.motifs.isEmpty()) {
                    ppr.mergeMotifs();
                    ppr.getBestString(); // calculate best_sob;
                    score = ppr.best_sob != null ? ppr.best_sob.getScore() : 0.0;
                    if (score >= 40 && (ppr.best_sob.hasAdjacentMotifs() || ppr.best_sob.hasDYW())) {
                        beads_writer.write(ppr.getID() + "\t" + ppr.length() + "\t" + ppr.getBestString() + "\t"
                                + ppr.guessSubclass() + "\t" + score);
                        beads_writer.newLine();
                        orf_writer.write(">" + ppr.id);
                        orf_writer.newLine();
                        orf_writer.write(ppr.sequence);
                        orf_writer.newLine();
                        ppr.writeMotifs(motif_writer);
                    }
                }
            }

            orf_writer.close();
            motif_writer.close();
            beads_writer.close();

        }

        catch (Exception ex) {
            if (PPRfinder.VERBOSE)
                ex.printStackTrace();
            System.err.println(ex.toString());
            System.exit(-1);
        }
    }
}
