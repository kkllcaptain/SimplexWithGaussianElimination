/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package diet_with_gaussian;

import java.io.*;
import java.util.*;

/**
 *
 * @author sguo46
 */


class Equation {

    Equation(double a[][], double b[]) {
        this.a = a;
        this.b = b;
    }

    double a[][];
    double b[];
}

class Position {

    Position(int column, int row) {
        this.column = column;
        this.row = row;
    }

    int column;
    int row;
}


public class Diet_with_Gaussian {

    /**
     * @param args the command line arguments
     */
     BufferedReader br;
    PrintWriter out;
    StringTokenizer st;
    boolean eof;

    int solveDietProblem(int n, int m, double A[][], double[] b, double[] c, double[] x) {
        Arrays.fill(x, 1);
        // =============== Step 1. pre-process ==============================
        // take matrix A and vector b, build list A and List b
        // add non-negtivity constraints
        // remove parallel constraints
        // detect incorrect (non-overlapping) conditions and return -1
        // modify Alist and bList as they are objects
        // ==================================================================
        List<List> AList = new ArrayList<>();
        List<Double> bList = new ArrayList<>();
        boolean status = preProcess(A,b,AList,bList);
        if(!status) return -1;
        
        // =============== Step 2. Initialization Simplex Searching =========
        // Initialize empty lists that will store solutions, objective values
        // Make sets of subsets of combinations of equations that
        //                  each combination will be solved by Gaussian Elim
        // ==================================================================
                
        List<Double> OBJList = new ArrayList<>(); 
        List<Boolean> Valid = new ArrayList<>();
        List<double[]> SolnList = new ArrayList<>();
        
        List<List> sets = setTraverse(bList.size(),AList.get(0).size());
        

        // ============== Step 3. Solving Each Combination ==================
        // ==================================================================
        for(int index =0; index<sets.size(); index++){
        // -------------- Step 3.a. rebuild new Equations -------------------
        // tak list AList, and bList and one subset of sets
        // generate Equation that contains matrix a and vector b
        // ------------------------------------------------------------------   

        List<Integer> subset = sets.get(index);
        Equation myEQ = buildEquations(AList, bList, subset);
        
        
        
        
        
        
        
        }
        
        

        
        return 0;
    }
      
    static boolean preProcess(double[][] A, double[] b, List<List> AList, List<Double> bList){
        // ==================== Build Initial A,b Lists =====================
        // Matrix (Arrays) A and b can not be modified, transfer these to 
        // arrayLists first
        // ==================================================================
        for(int i =0; i<A.length; i++){
            List RowList = new ArrayList<Double>();
            for(int j=0; j<A[0].length; j++){
                RowList.add(A[i][j]);
            }
            AList.add(new ArrayList(RowList));
            bList.add(b[i]);
        }

//        System.out.println("initial");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
        
        
        
        // =============== add non-negativity constraints ===================
        // all x's >= 0   this should be written as    -x_i <= 0
        // ==================================================================
        int n_x = AList.get(0).size();
        for(int i=0; i<n_x; i++){
            List newRow = new ArrayList<Double>();
            for(int j=0; j<n_x; j++){
                double element = (j==i) ? -1.0: 0.0;
                newRow.add(element);
            }
            AList.add(new ArrayList(newRow));
            bList.add((double) 0);
        }
        
//        System.out.println("after adding non-negativity");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
        
        // ==================== Check for Parallel Constraints ==============
        // if A!/A       ==> it is fine, not parallel 
        // if A//A, b//b ==> remove one
        // if A//A, b!/b ==> case 1: one constraint is redundant and should be 
        //                           removed
        //                   case 2: there is no-over lapping of the two 
        //                           ==> report no solution
        //                   case 3: over lapping, that's fine, do nothing
        // ==================================================================
        for(int row =0; row<AList.size(); row++){for(int otherRow = 0; otherRow<AList.size(); otherRow++){if(row != otherRow){
            // ratios is the ratio between each pair of elements of rows in 
            // matrix A 
            //System.out.println("comparing row: "+row+" and other row: "+otherRow);
            List<Double> ratios = new ArrayList<>(); 
            double ARatio = (double) 0.0;
            double bRatio = (double) 0.0;
            
            for(int col =0; col<AList.get(0).size(); col++){
                double tempRatio = (double)AList.get(row).get(col)/(double)AList.get(otherRow).get(col);
                boolean isZero = Objects.equals((double)AList.get(row).get(col), 0.0);
                if (!(Double.isInfinite(tempRatio)||Double.isNaN(tempRatio))||!isZero){
                    ratios.add(tempRatio);
                    //System.out.println(tempRatio);
                }
            }
                bRatio = (double)bList.get(row)/(double)bList.get(otherRow);
            
            if(ratios.stream().distinct().limit(2).count() <= 1){ // A parallel
                //System.out.println("A parallel");
                ARatio = (double) ratios.get(0);
                if(Objects.equals(ARatio, bRatio)){
                    //System.out.println("b parallel");
                    AList.remove(otherRow);
                    bList.remove(otherRow);
                    //System.out.println("*******************removing row: "+otherRow);
                } else if(ARatio > 0){
                    //System.out.println("same direction");
                    if (bList.get(otherRow) > bList.get(row)) {
                        AList.remove(otherRow);
                        bList.remove(otherRow);
                        //System.out.println("*******************removing row: "+otherRow);
                    } else {
                        AList.remove(row);
                        bList.remove(row);
                        //System.out.println("*******************removing row: "+row);
                    }
                } else if ((((double)bList.get(row))+((double)bList.get(otherRow)))<0){
                    //System.out.println("not overlapping");
                    return false; // answer if the system of equations are fine
                }
            }
            
//        System.out.println("after processing");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        System.out.println();
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
            
            
        }}}
        
        return true;
    }
    
    static boolean checkSolution(double[] solution, Equation myNewEQ) {
        
        double[][] A = myNewEQ.a;
        double[] b = myNewEQ.b;
        
        for (int i = 0; i < A.length; i++) {
            double[] inequality = A[i];
//            System.out.println("solution: "+ Arrays.toString(solution));
//            System.out.println("inequality: "+ Arrays.toString(inequality));
//            System.out.println("b: "+ b[i]);
//            System.out.println("sum: "+ Arrays.stream(multiplyArrays(solution, inequality)).sum() );
            if (Arrays.stream(multiplyArrays(solution, inequality)).sum() > b[i])
                return false;
        }
        return true;
    }
    
    static double[] multiplyArrays(double[] a, double b[]) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] * b[i];
        }
        return result;
    }
    
    
    static List<List> setTraverse(int n, int m){
        List<List> Sets = new ArrayList<>();
        // n # of constraints/restrictions
        // m # of variables
        // m out of n constrains should solve to an point
        if(m == n){
            List<Integer> subset = new ArrayList<>();
            for (int iter = 0; iter<m; iter ++){
                subset.add(1);
            }
            Sets.add(subset);
            return Sets;
        } else {

            for (int i = 0; i < (1 << n); i++) {
                ArrayList<Integer> subset = new ArrayList<>();
                for (int j = 0; j < n; j++) {
                    if (((i >> j) & 1) == 1) {
                        subset.add(1);
                    } else subset.add(0);
                }
                
                int sum = subset.stream().mapToInt(Integer::intValue).sum();
                if (sum == m) Sets.add(subset);                
                
            }

            return Sets;
        }
    }
    
    static Equation buildEquations(List<List> AList, List<Double> bList, List<Integer> subset){
        List<ArrayList> AE = new ArrayList(); 
        List<Double> bE = new ArrayList(); 
        
        for(int index =0; index<AList.size(); index++){
            if((int)subset.get(index) == 1){
                AE.add(new ArrayList((AList.get(index))));
                bE.add(bList.get(index));
            }
        }

        double[][] Ap = new double[AE.size()][AE.get(0).size()];
        double[] bp = new double[AE.size()];
        
        for (int row = 0; row < AE.size(); row++){
            bp[row] = (double) bE.get(row);
            for (int column = 0; column < AE.get(0).size(); column ++){
                Ap[row][column] = (double) AE.get(row).get(column);
            }
            
        }
        
        Equation myEq = new Equation(Ap,bp);
        
        return myEq;
       
    }
    
    void solve() throws IOException {
        int n = nextInt();
        int m = nextInt();
        double[][] A = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = nextInt();
            }
        }
        double[] b = new double[n];
        for (int i = 0; i < n; i++) {
            b[i] = nextInt();
        }
        double[] c = new double[m];
        for (int i = 0; i < m; i++) {
            c[i] = nextInt();
        }
        double[] ansx = new double[m];
        int anst = solveDietProblem(n, m, A, b, c, ansx);
        if (anst == -1) {
            out.printf("No solution\n");
            return;
        }
        if (anst == 0) {
            out.printf("Bounded solution\n");
            for (int i = 0; i < m; i++) {
                out.printf("%.18f%c", ansx[i], i + 1 == m ? '\n' : ' ');
            }
            return;
        }
        if (anst == 1) {
            out.printf("Infinity\n");
            return;
        }
    }

    Diet_with_Gaussian () throws IOException {
        br = new BufferedReader(new InputStreamReader(System.in));
        out = new PrintWriter(System.out);
        solve();
        out.close();
    }

    public static void main(String[] args) throws IOException {
        new Diet_with_Gaussian();
    }

    String nextToken() {
        while (st == null || !st.hasMoreTokens()) {
            try {
                st = new StringTokenizer(br.readLine());
            } catch (Exception e) {
                eof = true;
                return null;
            }
        }
        return st.nextToken();
    }

    int nextInt() throws IOException {
        return Integer.parseInt(nextToken());
    }
    
}
