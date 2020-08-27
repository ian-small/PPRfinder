/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

/**
 *
 * @author ian
 */
public class Motif {

    static final int MAX_MOTIF_LENGTH = 60;
    static final int CANONICAL_P_LENGTH = 35;
    static final int CANONICAL_P1_LENGTH = 35;
    static final int CANONICAL_P2_LENGTH = 35;
    static final int CANONICAL_S1_LENGTH = 31;
    static final int CANONICAL_S2_LENGTH = 32;
    static final int CANONICAL_SS_LENGTH = 31;
    static final int CANONICAL_L1_LENGTH = 35;
    static final int CANONICAL_L2_LENGTH = 36;
    static final int CANONICAL_E1_LENGTH = 34;
    static final int CANONICAL_E2_LENGTH = 34;
    static final int CANONICAL_DYW_LENGTH = 136;
    static final int CANONICAL_KPAKA_LENGTH = 135;

    public int hmm_start;
    public int hmm_end;
    public int alignment_start;
    public int alignment_end;
    public int inferred_start;
    public int inferred_end;
    public double score;
    public String hmmer_alignment = null;
    public String type;
    public boolean hidden = false; // flag to decide whether or not to include motif (because it has a marginal
                                   // score)

    public Motif(String s) {
        // dummy constructor for blank Motif to act as a gap or end
        type = s;
    }

    static int getCanonicalLength(Motif m) {
        switch (m.type) {
            case "P":
                return CANONICAL_P_LENGTH;
            case "P1":
                return CANONICAL_P1_LENGTH;
            case "P2":
                return CANONICAL_P2_LENGTH;
            case "S1":
                return CANONICAL_S1_LENGTH;
            case "S2":
                return CANONICAL_S2_LENGTH;
            case "SS":
                return CANONICAL_SS_LENGTH;
            case "L1":
                return CANONICAL_L1_LENGTH;
            case "L2":
                return CANONICAL_L2_LENGTH;
            case "E1":
                return CANONICAL_E1_LENGTH;
            case "E2":
                return CANONICAL_E2_LENGTH;
            case "DYW":
                return CANONICAL_DYW_LENGTH;
            case "KPAKA":
                return CANONICAL_KPAKA_LENGTH;
            default:
                return -1;
        }
    }

    public Motif(String[] fields) {

        type = fields[3];
        score = Double.parseDouble(fields[13]);
        hmm_start = Integer.parseInt(fields[15]);
        hmm_end = Integer.parseInt(fields[16]);
        alignment_start = Integer.parseInt(fields[17]);
        alignment_end = Integer.parseInt(fields[18]);

    }

    public int getHMMstart() {
        return hmm_start;
    }

    public void setAlignment(String alignment) {
        // remove leading and trailing gaps
        int i = 0;

        while (alignment.charAt(i) == '-') {
            i++;
        }

        int j = alignment.length() - 1;

        while (alignment.charAt(j) == '-') {
            j--;
        }

        hmmer_alignment = alignment.substring(i, j + 1);
    }

    public boolean inferStart(int limit) {
        inferred_start = alignment_start - (hmm_start - 1);
        if (inferred_start < limit) {
            // System.out.println("setting motif start to limit for motif " + type);
            inferred_start = limit;
            return false;
        } else {
            return true;
        }
    }

    public void inferEnd(int limit) {
        inferred_end = limit;
    }

    public int get_motif_length() {
        return hmm_end + get_C_extension();
    }

    public double hmm_coverage() {
        // return ((double)(hmm_end-hmm_start+1))/getCanonicalLength(this); //total
        // coverage
        return ((double) (hmm_end)) / getCanonicalLength(this); // C-ter coverage
    }

    public int get_N_extension() {
        return alignment_start - inferred_start;
    }

    public int get_C_extension() {
        return inferred_end - alignment_end;
    }

    public String getType() {
        switch (type) {
            case "DYW":
                if (hmm_coverage() < 0.97) {
                    return "E+";
                } else {
                    return "DYW";
                }
            case "KPAKA":
                if (hmm_coverage() < 0.97) {
                    return "E+:KP";
                } else {
                    return "DYW:KP";
                }
            default:
                return type;
        }
    }

    public String toString() {
        return type + "\t" + String.valueOf(alignment_start) + "\t" + String.valueOf(alignment_end) + "\t"
                + String.valueOf(score);
    }

}
