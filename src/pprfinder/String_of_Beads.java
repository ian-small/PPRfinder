/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.util.ArrayList;

/**
 *
 * @author ian
 */
public class String_of_Beads {

    ArrayList<Motif> motifs;

    public String_of_Beads(Motif m) {

        motifs = new ArrayList<Motif>();
        motifs.add(m);

    }

    public String_of_Beads(ArrayList<Motif> string, Motif m) {

        motifs = string;
        motifs.add(m);

    }

    public boolean canAddMotif(Motif m) {

        if (m.hidden)
            return false;

        if (getLastMotif().alignment_end >= m.alignment_start) {
            return false;
        } else {
            return true;
        }
    }

    public void addMotif(Motif m) {
        motifs.add(m);
    }

    public ArrayList<Motif> getMotifs() {
        return new ArrayList<Motif>(motifs); // more typesafety
        // return (ArrayList<Motif>) motifs.clone();
    }

    public Motif getLastMotif() {
        return motifs.get(motifs.size() - 1);
    }

    public boolean hasAdjacentMotifs() {
        if (motifs.size() < 2)
            return false;
        Motif m1, m2;
        for (int i = 0; i < (motifs.size() - 2); i++) {
            m1 = motifs.get(i);
            m2 = motifs.get(i + 1);
            if (m2.inferred_start == m1.inferred_end + 1)
                return true;
        }
        return false;
    }

    public boolean hasDYW() {
        Motif m;
        for (int i = 0; i < motifs.size(); i++) {
            m = motifs.get(i);
            if (m.type.startsWith("DYW") || m.type.equals("KPAKA"))
                return true;
        }
        return false;
    }

    // public ArrayList<Motif> getSubstring() {
    // ArrayList<Motif> substring = new ArrayList<Motif>();

    // for (int i=0; i < (motifs.size()-1); i++) {
    // substring.add(motifs.get(i));
    // }

    // return substring;
    // }

    public double getScore() {
        double total = 0.0;
        Motif m;
        for (int i = 0; i < motifs.size(); i++) {
            m = (Motif) motifs.get(i);
            total += m.score;
        }
        return total;
    }

    Motif inferiorMotif(Motif m1, Motif m2) {
        if (m1.type.equals("DYW") && m2.type.equals("KPAKA")) {
            if (-0.45522599 * m1.score + 0.52263719 * m2.score - 5.28553775 > 0) { // SVM model
                return m1;
            }
            return m2;
        }
        if (m1.type.equals("KPAKA") && m2.type.equals("DYW")) {
            if (-0.45522599 * m2.score + 0.52263719 * m1.score - 5.28553775 > 0) {
                return m2;
            }
            return m1;
        }
        if (m1.score > m2.score) {
            return m2;
        } else {
            return m1;
        }
    }

    public void inferMotifs(PPR ppr) {
        int limit;

        Motif m1, m2;

        boolean overlap;

        do {
            overlap = false;
            // infer starts
            for (int i = motifs.size() - 1; i >= 0; i--) {

                m2 = motifs.get(i);

                if (i > 0) { // infer start with respect to preceding motif
                    m1 = motifs.get(i - 1);
                    limit = m1.alignment_end + 1;
                    if (!m2.inferStart(limit)) { // motifs overlap, so remove lowest scoring match
                        motifs.remove(inferiorMotif(m1, m2));
                        overlap = true;
                        break;
                    }
                } else {
                    limit = 1;
                }

                // System.out.println("inferring motif start: " + m2.type + "\t" +
                // String.valueOf(m2.inferred_start));
                m2.inferStart(limit);

                if (m2.inferred_start > m2.alignment_start) {
                    System.err.println("inferred start exceeds alignment start for " + ppr.getID() + " motif "
                            + String.valueOf(i + 1));
                }
            }
        } while (overlap);

        int canonical_length;

        // infer ends
        for (int i = 0; i < motifs.size(); i++) {

            boolean define_end_by_start_of_next_motif = false;

            limit = ppr.length();

            m1 = motifs.get(i);
            if (i + 1 < motifs.size()) {
                m2 = motifs.get(i + 1);
                if (m1.inferred_start + Motif.MAX_MOTIF_LENGTH > m2.inferred_start) {
                    define_end_by_start_of_next_motif = true;
                    limit = m2.inferred_start - 1;
                    // System.out.println("inferring motif end: " + m1.type + " " +
                    // String.valueOf(m1.inferred_end) + " with respect to motif " + m2.type + "
                    // starting at " + String.valueOf(m2.inferred_start) );
                }
            }

            if (!define_end_by_start_of_next_motif) {

                canonical_length = Motif.getCanonicalLength(m1);

                switch (m1.getType()) {
                    case "P":
                    case "P1":
                    case "P2":
                    case "S1":
                    case "S2":
                    case "SS":
                    case "L1":
                    case "L2":
                    case "E1":
                    case "E2":
                        limit = Math.min(limit, m1.inferred_start + canonical_length - 1);
                        break;
                    case "E+":
                    case "DYW":
                    case "KPAKA":
                        // leave limit at end of sequence
                        break;
                    default:
                        break;
                }

            }

            m1.inferEnd(limit);

        }
    }

    public boolean isValid() {

        Motif m1, m2, m3;
        Boolean b;

        for (int i = 0; i < motifs.size() - 1; i++) {

            m1 = motifs.get(i);
            m2 = (Motif) motifs.get(i + 1);

            String combo = m1.type + "." + m2.type;
            // System.out.println(combo);

            b = PPRfinder.Gapless_Motif_Combinations.get(combo);

            if (b == null) {
                if (PPRfinder.VERBOSE)
                    System.out.println("No Gapless motif combination for: " + combo);
                return false;
            }

            if (!b) {
                return false;
            }
        }

        Motif gap = new Motif("-");
        Motif end = new Motif("*");

        for (int i = 0; i < motifs.size(); i++) {

            if (i == 0) {
                m1 = end;
            } else {
                m1 = motifs.get(i - 1);
            }
            m2 = (Motif) motifs.get(i);

            if (i + 1 == motifs.size()) {
                m3 = end;
            } else {
                m3 = (Motif) motifs.get(i + 1);
            }

            if (m1 != gap && m1.inferred_end + 30 < m2.inferred_start) {
                m3 = m2;
                m2 = gap;
                if (m1 != gap && m1.inferred_end + 60 < m3.inferred_start) {
                    m3 = gap; // double gap
                }
            }

            if (m2 != gap && m3 != gap && m2.inferred_end + 8 < m3.inferred_start) {
                m3 = gap;
            }

            String combo = m1.type + "." + m2.type + "." + m3.type;

            if (!PPRfinder.Motif_Combinations.get(combo)) {
                return false;
            }
        }

        return true;
    }

    public String beads(PPR ppr) {
        StringBuilder buff = new StringBuilder();
        Motif m = null;
        int p = 1;
        int gap;
        for (int i = 0; i < motifs.size(); i++) {
            m = (Motif) motifs.get(i);
            if (m == null) {
                System.err.println(i + " null motif");
                continue;
            }
            gap = m.inferred_start - p;
            if (gap > 0) {
                buff.append(String.valueOf(gap)).append("-");
            }
            p = m.inferred_end + 1;
            buff.append(m.getType());
            if (i < motifs.size() - 1) {
                buff.append("-");
            }
        }
        gap = ppr.sequence.length() - m.inferred_end;
        if (gap > 0) {
            buff.append("-").append(String.valueOf(gap));
        }
        return buff.toString();
    }

    public boolean hasMotifType(String mt) {
        for (Motif m : motifs) {
            if (m.getType().equals(mt)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String toString() {
        StringBuilder string = new StringBuilder();
        for (Motif m : motifs) {
            string.append(m.getType());
        }
        return string.toString();
    }
}
