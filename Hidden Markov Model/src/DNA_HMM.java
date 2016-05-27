
/**
 * Ainshames University, Faculty of Engineering
 * Computer Engineering and Software Systems department
 * @author Aly Osama Aly
 * @email alyosamah@gmail.com
 */
import com.mxgraph.layout.mxCompactTreeLayout;
import com.mxgraph.layout.mxIGraphLayout;
import com.mxgraph.swing.mxGraphComponent;

import java.awt.HeadlessException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import org.jgrapht.ListenableGraph;
import org.jgrapht.ext.JGraphXAdapter;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.ListenableDirectedWeightedGraph;


public class DNA_HMM extends javax.swing.JFrame {

    public enum STATE_TYPE {
        INSERT, DELETE, MATCH, BEGIN, END
    };

    public enum DNA {
        A, C, G, T, GAP
    };

    final int M_M = 0;
    final int M_D = 1;
    final int M_I = 2;

    final int I_M = 0;
    final int I_I = 1;

    final int D_M = 0;
    final int D_D = 1;

    private static double[][] emissionM;
    private static double[][] emissionI;
    private static double[][] transM;
    private static double[][] transI;
    private static double[][] transD;
    private  static int[][] prob;

    private int nCons = 0;

    /**
     * Creates new form DNA_HMM
     */
    public DNA_HMM() {
        initComponents();
    }

    /**
     * Build HMM Model From MSA Generates Emission and Transition Matrices
     *
     * @param MSA Multiple Sequence Alignment
     * @param sizeOfData number of sequences
     * @param sizeOfSeq Size of Characters in any alignment
     */
    public void buildHMM(char[][] MSA, int sizeOfData, int sizeOfSeq) {
        ////Step 2: Check and Find Conservative columns
        int[] conservateCol = new int[sizeOfSeq];
        nCons=0;
        for (int i = 0; i < sizeOfSeq; i++) {
            for (int j = 0; j < 4; j++) {
                if (prob[j][i] > sizeOfData / 2) {
                    nCons += 1;
                    conservateCol[i] = nCons;

                }
            }
        }

    
        ///////////////// Step 3: Filling Emission and Transition Matrices \\\\\\\\\\\\\\\\\\\\
        //// Intiallizations 
        emissionM = create2Darray(4, nCons + 1);
        emissionI = create2Darray(4, nCons + 1);
        transM = create2Darray(3, nCons + 1);
        transI = create2Darray(2, nCons + 1);
        transD = create2Darray(2, nCons + 1);

        int lastConservative = 0;

        //// Filling Emission Matrices (Counting)
        for (int i = 0; i < sizeOfSeq; i++) {
            // Check if conservative column
            if (conservateCol[i] != 0) {
                lastConservative = conservateCol[i];
                //Fill Emission of Matching state
                for (int k = 0; k < 4; k++) {
                    emissionM[k][conservateCol[i]] = prob[k][i];
                }
            } else {
                //Fill Emission of Insertion state
                for (int k = 0; k < 4; k++) {
                    emissionI[k][lastConservative] += prob[k][i];
                }
            }
        }

        //// Filling Transtion Matrices (Counting)
        for (int i = 0; i < sizeOfData; i++) {
            lastConservative = 0;
            STATE_TYPE state = STATE_TYPE.BEGIN;
            for (int j = 0; j < sizeOfSeq; j++) {
                switch (state) {
                    case BEGIN:
                        if (conservateCol[j] != 0) {
                            lastConservative = conservateCol[j];
                            if (MSA[i][j] == '-') {
                                transM[M_D][0] += 1;
                                state = STATE_TYPE.DELETE;
                            } else {
                                transM[M_M][0] += 1;
                                state = STATE_TYPE.MATCH;
                            }
                        } else if (MSA[i][j] == '-') {
                            state = STATE_TYPE.MATCH;
                        } else {
                            transM[M_I][0] += 1;
                            state = STATE_TYPE.INSERT;
                        }
                        j -= 1;
                        break;
                    case MATCH:
                        if (conservateCol[j + 1] != 0) {
                            lastConservative = conservateCol[j + 1];
                            if (MSA[i][j + 1] != '-') {
                                transM[M_M][conservateCol[j + 1] - 1] += 1;
                                state = STATE_TYPE.MATCH;
                            } else {
                                transM[M_D][conservateCol[j + 1] - 1] += 1;
                                state = STATE_TYPE.DELETE;
                            }
                        } else if (MSA[i][j + 1] != '-') {
                            transM[M_I][lastConservative] += 1;
                            state = STATE_TYPE.INSERT;
                        }

                        if (j == sizeOfSeq - 2) {
                            state = STATE_TYPE.END;
                        }
                        break;

                    case INSERT:
                        if (conservateCol[j + 1] != 0) {
                            lastConservative = conservateCol[j + 1];
                            if (MSA[i][j + 1] != '-') {
                                transI[I_M][conservateCol[j + 1] - 1] += 1;
                                state = STATE_TYPE.MATCH;

                            } else if (MSA[i][j + 1] != '-') {
                                transI[I_I][lastConservative] += 1;
                                state = STATE_TYPE.INSERT;
                            }
                        } else if (MSA[i][j + 1] != '-') {
                            transI[I_I][lastConservative] += 1;
                            state = STATE_TYPE.INSERT;
                        }
                        if (j == sizeOfSeq - 2) {
                            state = STATE_TYPE.END;
                        }
                        break;

                    case DELETE:
                        if (conservateCol[j + 1] != 0) {
                            lastConservative = conservateCol[j + 1];
                            if (MSA[i][j + 1] != '-') {
                                transD[D_M][conservateCol[j + 1] - 1] += 1;
                                state = STATE_TYPE.MATCH;

                            } else if (MSA[i][j + 1] != '-') {
                                transD[D_D][lastConservative] += 1;
                                state = STATE_TYPE.DELETE;
                            }
                        }

                        if (j == sizeOfSeq - 2) {
                            state = STATE_TYPE.END;
                        }
                        break;

                    case END:
                        if (conservateCol[j] != 0) {
                            if (MSA[i][j] == '-') {
                                transM[M_D][conservateCol[j]] += 1;
                                state = STATE_TYPE.DELETE;
                            } else {
                                transM[M_M][conservateCol[j]] += 1;
                                state = STATE_TYPE.MATCH;
                            }
                        } else if (MSA[i][j] != '-') {
                            transM[M_I][conservateCol[j]] += 1;
                            state = STATE_TYPE.INSERT;
                        }
                        break;
                    default:
                        break;
                }

            }
        }

        //// (Computing Probabilities ) Filling Emission and Transition Matrices 
        emissionM = computeProbaility(emissionM);
        emissionI = computeProbaility(emissionI);
        transM = computeProbaility(transM);
        transI = computeProbaility(transI);
        transD = computeProbaility(transD);

        ////////////////// PRINT OUTPUT /////////////////////////
        printOutput();
    }

    /***
     * 
     * Drawing HMM Model Graph to the Screen
     * @return 
     */
    public ListenableGraph<String, MyEdge> buildGraph() {
        
        ListenableDirectedWeightedGraph<String, MyEdge> g
                = new ListenableDirectedWeightedGraph<>(MyEdge.class);

        
        String x0 = "start";
        String xend = "end";
        String[] M = new String[nCons];
        String[] D = new String[nCons];
        String[] I = new String[nCons + 1];

        g.addVertex(x0);
        g.addVertex(xend);
        MyEdge e;

        String[] state_str = {"M", "D", "I"};        
        /// Drawing Matching States
        for (int i = 0; i < nCons; i++) {
            if (transM[M_M][i] != 0) {
                M[i] = state_str[M_M] + String.valueOf(i + 1);
                g.addVertex(M[i]);
                if (i == 0) {
                    e = g.addEdge(x0, M[i]);
                } else {
                    e = g.addEdge(M[i - 1], M[i]);
                }
                g.setEdgeWeight(e, transM[M_M][i]);
                if (i == nCons - 1) {
                    e = g.addEdge(M[i], xend);
                    g.setEdgeWeight(e, transM[M_M][nCons]);
                }
            }
        }

        /// Drawing Insertion States
        for (int i = 0; i <= nCons; i++) {
            if (transM[M_I][i] != 0) {
                I[i] = state_str[M_I] + String.valueOf(i);
                g.addVertex(I[i]);

                if (i == 0) {
                    e = g.addEdge(x0, I[i]);
                } else {
                    e = g.addEdge(M[i - 1], I[i]);
                }
                g.setEdgeWeight(e, transM[M_I][i]);
                if (i == nCons) {
                    e = g.addEdge(I[i], xend);
                    g.setEdgeWeight(e, transM[M_I][nCons]);
                }

                if (transI[I_M][i] != 0) {
                    e = g.addEdge(I[i], M[i]);
                    g.setEdgeWeight(e, transI[I_M][i]);
                }
                if (transI[I_I][i] != 0) {
                    e = g.addEdge(I[i], I[i]);
                    g.setEdgeWeight(e, transI[I_I][i]);
                }
            }
        }

        /// Drawing Deletion States
        for (int i = 0; i <= nCons; i++) {
            if (transM[M_D][i] != 0) {
                D[i] = state_str[M_D] + String.valueOf(i + 1);
                g.addVertex(D[i]);
                if (i == 0) {
                    e = g.addEdge(x0, D[i]);
                } else {
                    e = g.addEdge(M[i - 1], D[i]);
                }
                g.setEdgeWeight(e, transM[M_D][i]);
                if (i == nCons - 1) {
                    e = g.addEdge(D[i], xend);
                    g.setEdgeWeight(e, transM[M_D][nCons]);
                }

                if (transD[D_M][i + 1] != 0) {
                    e = g.addEdge(D[i], M[i]);
                    g.setEdgeWeight(e, transD[D_M][i + 1]);
                }
                if (transD[D_D][i + 1] != 0) {
                    e = g.addEdge(D[i], D[i]);
                    g.setEdgeWeight(e, transD[D_D][i + 1]);
                }
            }
        }

        return g;
    }

    public static class MyEdge extends DefaultWeightedEdge {

        @Override
        public String toString() {
            return String.valueOf(getWeight());
        }
    }

    /**
     * 
     * Counts every Nucleotide in the multiple sequence alignment
     * @param MSA
     * @return 
     */
    public boolean countN(char[][] MSA){
        int sizeOfSeq=MSA[0].length;
        int sizeOfData=MSA.length;
       
        ////Step 1:   Count Number of Letters in each column of MSA
        prob = create2DarrayInt(DNA.values().length, sizeOfSeq);
        for (int i = 0; i < sizeOfSeq; i++) {
            for (int j = 0; j < sizeOfData; j++) {
                switch (MSA[j][i]) {
                    case 'A':
                        prob[DNA.A.ordinal()][i] += 1;
                        break;
                    case 'C':
                        prob[DNA.C.ordinal()][i] += 1;
                        break;
                    case 'G':
                        prob[DNA.G.ordinal()][i] += 1;
                        break;
                    case 'T':
                        prob[DNA.T.ordinal()][i] += 1;
                        break;
                    case '-':
                    case '_':
                        prob[DNA.GAP.ordinal()][i] += 1;
                        break;
                    default:
                        JOptionPane.showMessageDialog(null, "Invalid Inputs not A/C/G/T/-", "Input Error", JOptionPane.WARNING_MESSAGE);
                        return false;
                }

            }
        }
        return true;
    }
    
    /**
     * 
     * @param array array that we want to print
     * @param message array name
     * @param rowsNames representation of every column
     */
    public void print2DArray(double[][] array,String message,String[] rowsNames){
        
        int x=array.length;
        int y=array[0].length;
        outputTxtArea.append("\n"+message+" Table:\n");
        for(int i=0;i<y;i++)outputTxtArea.append("\t"+i+" ");
        outputTxtArea.append("\n");
        for (int i = 0; i < x; i++) {
            outputTxtArea.append(rowsNames[i] + "\t");
            for (int j = 0; j < y; j++) {
                outputTxtArea.append(array[i][j] + "\t");
            }
            outputTxtArea.append("\n");
        }
    }
    /**
     * Print the transition and emission matrices to screen
     */
    public void printOutput() {
 
        print2DArray(emissionM,"Emission Match",new String[]{"A", "C", "G", "T"});
        print2DArray(emissionI,"Emission Insertion",new String[]{"A", "C", "G", "T"});
        print2DArray(transM,"Transition MATCH",new String[]{"M_M", "M_D", "M_I"});
        print2DArray(transI,"Transition INSERTION",new String[]{"I_M", "I_I"});
        print2DArray(transD,"Transition DELETION",new String[]{"D_M", "D_D"});
         /*
        // Print 
        for (int i = 0; i < sizeOfSeq; i++) {
            if (conservateCol[i] != 0) {
                outputTxtArea.append(conservateCol[i] + "-" + i + " ");
            }
        }
        outputTxtArea.append("\n\n");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < sizeOfSeq; j++) {
                outputTxtArea.append(prob[i][j] + " ");
            }
            outputTxtArea.append("\n");
        }
         */
    }

    /////////////////////// HELPER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\
    double d2p(double x) {
        return (double) Math.round(x * 100d) / 100d;
    }

    double[][] create2Darray(int n, int m) {
        double[][] array = new double[n][];
        for (int i = 0; i < n; i++) {
            array[i] = new double[m];
        }
        return array;
    }

    int[][] create2DarrayInt(int n, int m) {
        int[][] array = new int[n][];
        for (int i = 0; i < n; i++) {
            array[i] = new int[m];
        }
        return array;
    }

    double[][] computeProbaility(double[][] array) {
        int x = array.length;
        int y = array[0].length;
        for (int j = 0; j < y; j++) {
            double sum = 0;
            for (int i = 0; i < x; i++) {
                sum += array[i][j];
            }
            for (int i = 0; i < x; i++) {
                if (sum != 0) {
                    array[i][j] /= sum;
                    array[i][j] = d2p(array[i][j]);
                }
            }
        }
        return array;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        msaTxtArea = new javax.swing.JTextArea();
        loadBtn = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        buildBtn = new javax.swing.JButton();
        jScrollPane2 = new javax.swing.JScrollPane();
        outputTxtArea = new javax.swing.JTextArea();
        jLabel2 = new javax.swing.JLabel();
        saveBtn = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("DNA MSA");

        msaTxtArea.setColumns(20);
        msaTxtArea.setFont(new java.awt.Font("Monospaced", 0, 16)); // NOI18N
        msaTxtArea.setRows(5);
        jScrollPane1.setViewportView(msaTxtArea);

        loadBtn.setText("Load from File");
        loadBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadBtnActionPerformed(evt);
            }
        });

        jLabel1.setText("(Input) Multiple DNA Sequence Alignment : ");

        buildBtn.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        buildBtn.setText("Build HMM");
        buildBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buildBtnActionPerformed(evt);
            }
        });

        outputTxtArea.setEditable(false);
        outputTxtArea.setColumns(20);
        outputTxtArea.setFont(new java.awt.Font("Monospaced", 0, 14)); // NOI18N
        outputTxtArea.setRows(5);
        jScrollPane2.setViewportView(outputTxtArea);

        jLabel2.setText("(Output) Transition and Emission Tables: ");

        saveBtn.setText("Save to File");
        saveBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveBtnActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 334, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(buildBtn))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(loadBtn)))
                .addGap(8, 8, 8)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addGap(18, 18, 18)
                        .addComponent(saveBtn)
                        .addContainerGap())
                    .addComponent(jScrollPane2)))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(loadBtn)
                    .addComponent(jLabel2)
                    .addComponent(saveBtn))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 229, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 229, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(88, 88, 88)
                        .addComponent(buildBtn, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void loadBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadBtnActionPerformed
        fileChooser = new JFileChooser();
        int result = fileChooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            try {
                FileReader reader = new FileReader(file);
                try (BufferedReader br = new BufferedReader(reader)) {
                    msaTxtArea.read(br, null);
                }
                msaTxtArea.requestFocus();
            } catch (Exception e) {
            }
        }
    }//GEN-LAST:event_loadBtnActionPerformed

    private void saveBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveBtnActionPerformed
        fileChooser = new JFileChooser();
        int result = fileChooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            try {
                FileWriter writer = new FileWriter(file + ".txt");
                try (BufferedWriter bw = new BufferedWriter(writer)) {
                    outputTxtArea.write(bw);
                }
                //outputTxtArea.setText("");
                outputTxtArea.requestFocus();
                JOptionPane.showMessageDialog(null, "Saving the output done", "InfoBox: " + "Save Done", JOptionPane.INFORMATION_MESSAGE);

            } catch (IOException | HeadlessException e) {
            }
        }
    }//GEN-LAST:event_saveBtnActionPerformed

    private void buildBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buildBtnActionPerformed
        String[] inputLines = msaTxtArea.getText().split("\\n");
        
        //Check for input
        if (inputLines.length < 2) {
            JOptionPane.showMessageDialog(null, "Please write or load the input", "Input Error", JOptionPane.WARNING_MESSAGE);
        } else {

            //Check for equals alignment in input
            int sizeOfData = inputLines.length;
            char[][] MSA = new char[sizeOfData][];
            MSA[0] = inputLines[0].toUpperCase().toCharArray();
            int sizeOfSeq = MSA[0].length;
            for (int i = 1; i < inputLines.length; i++) {
                MSA[i] = inputLines[i].toUpperCase().toCharArray();
                if (MSA[i].length != sizeOfSeq) {
                    JOptionPane.showMessageDialog(null, "Incorrect Sequence Alignment (inequal lengths)", "Input Error", JOptionPane.WARNING_MESSAGE);
                    return;
                }
            }
            //Check for letters in input
            if(countN(MSA)){
                outputTxtArea.setText("");
                buildHMM(MSA, sizeOfData, sizeOfSeq);
                createAndShowGraph();
            }
        }


    }//GEN-LAST:event_buildBtnActionPerformed

    /**
     * 
     * Create the Frame of HMM Graph Model
    */
    private void createAndShowGraph() {
        JFrame frame = new JFrame("HMM Graph Model");
        ListenableGraph<String, MyEdge> g = buildGraph();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JGraphXAdapter<String, MyEdge> graphAdapter
                = new JGraphXAdapter<>(g);

        mxIGraphLayout layout = new mxCompactTreeLayout(graphAdapter);
        layout.execute(graphAdapter.getDefaultParent());

        frame.add(new mxGraphComponent(graphAdapter));

        frame.pack();
        frame.setLocationByPlatform(true);
        frame.setVisible(true);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(DNA_HMM.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }

        java.awt.EventQueue.invokeLater(() -> {
            new DNA_HMM().setVisible(true);
        });
    }

    private JFileChooser fileChooser;

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton buildBtn;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JButton loadBtn;
    private javax.swing.JTextArea msaTxtArea;
    private javax.swing.JTextArea outputTxtArea;
    private javax.swing.JButton saveBtn;
    // End of variables declaration//GEN-END:variables
}
