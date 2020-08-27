/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
// import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

class StartComparator implements Comparator<Motif> {
    @Override
    public int compare(Motif m1, Motif m2) {
        if (m1.alignment_start == m2.alignment_start) {
            return Double.compare(m2.score, m1.score);
        }
        return m1.alignment_start - m2.alignment_start;
    }
}

class TypeComparator implements Comparator<Motif> {
    @Override
    public int compare(Motif m1, Motif m2) {
        if (m1.type.equals(m2.type)) {
            return Integer.compare(m1.alignment_start, m2.alignment_start);
        } else
            return m1.type.compareTo(m2.type);
    }
}

class ScoreComparator implements Comparator<String_of_Beads> {
    @Override
    public int compare(String_of_Beads sob1, String_of_Beads sob2) {
        return Double.compare(sob1.getScore(), sob2.getScore());
    }
}

/**
 *
 * @author ian
 */
public class PPR {
    String id;
    String sequence;
    ArrayList<Motif> motifs;
    String_of_Beads best_sob;
    String_of_Beads best_valid_sob;

    public PPR(String name) {
        id = name;
        motifs = new ArrayList<Motif>();
    }

    public String getID() {
        return id;
    }

    public void setSequence(String data) {
        if (data.endsWith("*")) {
            sequence = data.substring(0, data.length() - 1); // removes trailing *
        } else {
            sequence = data.substring(0, data.length());
        }
    }

    public int length() {
        return sequence.length();
    }

    public boolean hasValidSob() {
        if (best_valid_sob == null) {
            return false;
        }
        return true;
    }

    public void addMotif(Motif m) {

        // bitscore filter
        if (m.type.equals("DYW") && m.score < 30)
            m.hidden = true;
        else if (m.type.equals("KPAKA") && m.score < 30)
            m.hidden = true;
        else if (m.type.equals("SS") && m.score < 10)
            m.hidden = true;
        else if (m.score < 0)
            m.hidden = true;

        if (!m.hidden)
            motifs.add(m);

    }

    public void addMotifAlignment(String motif_data, String motif_type) {
        String[] data = motif_data.split(" +");
        String[] coords = data[0].split("-");
        int coord_start = Integer.parseInt(coords[0]);
        int coord_end = Integer.parseInt(coords[1]);
        for (Motif m : motifs) {
            if (m.alignment_start == coord_start && m.alignment_end == coord_end && m.getType().equals(motif_type)) {
                m.setAlignment(data[1].replaceAll("-", ""));
                return;
            }
        }

        System.err.println("Could not find motif " + motif_data + " in " + id);
    }

    public void mergeMotifs() {

        if (motifs.size() < 2)
            return;

        Collections.sort(motifs, new TypeComparator());

        List<Motif> dead_motifs = new ArrayList<Motif>();
        List<Motif> live_motifs = new ArrayList<Motif>();

        Motif m1, m2;
        int m1_coverage, m2_coverage;
        double hmm_coverage, alignment_coverage;
        for (int i = 1; i < motifs.size(); i++) {
            m1 = motifs.get(i - 1);
            m2 = motifs.get(i);

            // check if merge is plausible
            if (!m1.type.equals(m2.type))
                continue;
            m1_coverage = m1.hmm_end - m1.hmm_start + 1;
            m2_coverage = m2.hmm_end - m2.hmm_start + 1;
            if (m1_coverage >= Motif.getCanonicalLength(m1))
                continue;
            if (m2_coverage >= Motif.getCanonicalLength(m2))
                continue;
            m2_coverage = m2.hmm_end - m2.hmm_start;
            hmm_coverage = m1_coverage + m2_coverage;
            alignment_coverage = m2.alignment_end - m1.alignment_start;
            if ((hmm_coverage / Motif.getCanonicalLength(m2) > 1.5)
                    || (alignment_coverage / Motif.getCanonicalLength(m2) > 1.5))
                continue;

            // merge
            if (PPRfinder.VERBOSE)
                System.out.println("merging " + id + " " + m2.type);
            Motif new_motif = new Motif(m2.type);
            new_motif.hmm_start = m1.hmm_start;
            new_motif.hmm_end = m2.hmm_end;
            new_motif.alignment_start = m1.alignment_start;
            new_motif.alignment_end = m2.alignment_end;
            new_motif.score = m1.score + m2.score;
            live_motifs.add(new_motif);
            dead_motifs.add(m1);
            dead_motifs.add(m2);
        }
        motifs.removeAll(dead_motifs);
        for (Motif m : live_motifs) {
            addMotif(m);
            // System.out.println("added " + id + " " + m.type);
        }

    }

    public String getBestString() {

        Collections.sort(motifs, new StartComparator());
        ArrayList<String_of_Beads> sobs = new ArrayList<String_of_Beads>();
        ArrayList<String_of_Beads> new_sobs = new ArrayList<String_of_Beads>();
        ArrayList<String_of_Beads> prunings = new ArrayList<String_of_Beads>();

        String_of_Beads new_sob;
        String_of_Beads extended_sob;

        for (Motif m : motifs) {
            new_sob = new String_of_Beads(m);
            new_sob.inferMotifs(this);
            if (new_sob.isValid()) {
                new_sobs.add(new_sob);
            }
            for (String_of_Beads sob : sobs) {
                if (sob.canAddMotif(m)) {
                    extended_sob = new String_of_Beads(sob.getMotifs(), m);
                    extended_sob.inferMotifs(this);
                    if (extended_sob.isValid()) {
                        new_sobs.add(extended_sob);
                    }
                }
            }
            sobs.addAll(new_sobs);
            new_sobs.clear();

            Collections.sort(sobs, new ScoreComparator());

            // prune sobs to avoid running out of memory...
            if (sobs.size() > 31) {
                String_of_Beads s;
                for (int n = 0; n < sobs.size() - 1; n++) {
                    s = sobs.get(n);
                    if (n < sobs.size() / 4) {
                        prunings.add(s);
                    }
                }
                sobs.removeAll(prunings);
                prunings.clear();
            }
            // System.out.println(id + " sobs contains " + sobs.size() + " members ");
        }

        /*
         * double best_score = -100.0; double best_valid_score = -100.0;
         * 
         * for (String_of_Beads sob : sobs) {
         * 
         * if (sob.getScore() > best_score) { best_sob = sob; best_score =
         * best_sob.getScore(); } if (sob.isValid() && sob.getScore() >
         * best_valid_score) { best_valid_sob = sob; best_valid_score =
         * best_valid_sob.getScore(); } }
         * 
         * if (best_valid_sob != null) { best_sob = best_valid_sob; }
         * best_sob.inferMotifs(this);
         */
        if (sobs.size() > 0) {
            best_sob = sobs.get(sobs.size() - 1);
            best_sob.inferMotifs(this);

            return best_sob.beads(this);
        }
        return "";
    }

    public int get_max_motif_length(String motif_type) {

        int max_motif_length = 0;

        for (Motif m : motifs) {
            if (m.hmmer_alignment != null && m.inferred_start <= m.alignment_start && m.getType().equals(motif_type)) {
                max_motif_length = Math.max(max_motif_length, m.get_motif_length());
            }
        }

        return max_motif_length;
    }

    public String guessSubclass() {

        if (best_sob.hasMotifType("P")) {
            return "P";
        }

        return "PLS";
    }

    /**
     *
     * @param motif_writer
     * @throws IOException
     */
    public void writeMotifs(BufferedWriter motif_writer) throws IOException {

        String_of_Beads sob = best_sob;

        try {

            for (Motif m : sob.getMotifs()) {
                motif_writer.write(id);
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.inferred_start));
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.inferred_end));
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.score));
                motif_writer.write("\t");
                motif_writer.write(sequence.substring(m.inferred_start - 1, m.inferred_end));
                motif_writer.write("\t");
                motif_writer.write(sequence.charAt(m.inferred_start));
                motif_writer.write("\t");
                motif_writer.write(sequence.charAt(m.inferred_start + 3));
                motif_writer.write("\t");
                motif_writer.write(sequence.charAt(m.inferred_end - 1));
                motif_writer.write("\t");
                motif_writer.write(m.getType());
                motif_writer.newLine();
            }

        } catch (StringIndexOutOfBoundsException ex) {
            System.err.println("StringIndexOutOfBoundsException for " + id);
        }

    }

    public void writeAlignments(BufferedWriter motif_writer, String motif_type, int max_motif_length)
            throws IOException {
        int hmmer_postgap_extension;

        String_of_Beads sob = best_sob;

        for (Motif m : sob.getMotifs()) {
            try {
                if (m.hmmer_alignment != null && m.inferred_start <= m.alignment_start
                        && m.getType().equals(motif_type)) {
                    motif_writer.write(">");
                    motif_writer.write(id);
                    motif_writer.write("/");
                    motif_writer.write(String.valueOf(m.inferred_start));
                    motif_writer.write("-");
                    motif_writer.write(String.valueOf(m.inferred_end));
                    motif_writer.newLine();

                    motif_writer.write(sequence.substring(m.inferred_start - 1, m.alignment_start - 1));

                    // write N-ter gaps (if any...)
                    for (int n = m.get_N_extension() + 1; n < m.hmm_start; n++) {
                        motif_writer.write(".");
                    }

                    // write hmmer alignment until 2nd loop gap site at position 29
                    hmmer_postgap_extension = m.hmm_end - 29;
                    if (hmmer_postgap_extension > 0) {
                        motif_writer.write(m.hmmer_alignment.substring(0,
                                m.hmmer_alignment.length() - hmmer_postgap_extension - 1));
                    } else {
                        motif_writer.write(m.hmmer_alignment);
                    }

                    for (int n = m.get_C_extension() + m.hmm_end; n < max_motif_length; n++) {
                        motif_writer.write(".");
                    }

                    if (hmmer_postgap_extension > 0) {
                        motif_writer.write(
                                m.hmmer_alignment.substring(m.hmmer_alignment.length() - hmmer_postgap_extension - 1));
                    }
                    motif_writer.write(sequence.substring(m.alignment_end, m.inferred_end));
                    motif_writer.newLine();
                }
            } catch (java.lang.StringIndexOutOfBoundsException e) {
                // ignore
            }
        }
    }

    // public double FFT(int wavelength, int phase) {

    // // values based on 100000000 random trials with Blosum62.random_associations
    // double random_mean = -0.9403338500000059;
    // double random_sd = 2.0946184866031654;

    // ArrayList<String>[] transforms = new ArrayList[wavelength];

    // for (int i = 0; i < transforms.length; i++) {
    // transforms[i] = new ArrayList<String>();
    // }

    // int index = phase;

    // while (index < sequence.length()) {
    // transforms[(index - phase) % wavelength].add(sequence.substring(index,
    // ++index));
    // }

    // double score = 0.0;
    // int count = 0;
    // for (ArrayList<String> transform : transforms) {
    // for (int i = 0; i < transform.size() - 1; i++) {
    // for (int j = i + 1; j < transform.size(); j++) {
    // score += Blosum62.score(transform.get(i), transform.get(j));
    // count++;
    // }
    // }
    // }

    // score /= count;

    // // convert to standard score
    // double z = (score - random_mean) / (random_sd / Math.sqrt(count));

    // return z;
    // }

}
