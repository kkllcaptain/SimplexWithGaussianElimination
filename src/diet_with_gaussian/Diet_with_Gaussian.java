/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//package diet_with_gaussian;

import java.io.*;
import java.util.*;

/**
 *
 * @author sguo46
 */


class Equation {
    List<ArrayList> a = new ArrayList<>();
    List<Double> b = new ArrayList<>();
    
    Equation(List<ArrayList> a, List<Double> b) {
        this.a = a;
        this.b = b;
    }


}

class Position {
    double value;
    int column;
    int row;

    Position(int column, int row) {
        this.column = column;
        this.row = row;
        this.value = 0.0;
    }
    
    Position(int column, int row, double value) {
        this.column = column;
        this.row = row;
        this.value = value;
    }
    

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
        // modify Alist and bList as they are objects
        // ==================================================================
        List<List> AList = new ArrayList<>();
        List<Double> bList = new ArrayList<>();
        preProcess(A,b,AList,bList);
        boolean noSoln = true;
        double maxOBJ = 0.0;
        int maxIndex = -1;
        
        // =============== Step 2. Initialization Simplex Searching =========
        // Initialize empty lists that will store solutions, objective values
        // Make sets of subsets of combinations of equations that
        //                  each combination will be solved by Gaussian Elim
        // ==================================================================
                
        List<Double> OBJList = new ArrayList<>(); 
        List<Boolean> Valid = new ArrayList<>();
        List<double[]> SolnList = new ArrayList<>();
        
        List<List> sets = setTraverse(bList.size(),AList.get(0).size());
        
        //System.out.println(Arrays.deepToString(sets.toArray()));

        // ============== Step 3. Solving Each Combination ==================
        // tak list AList, and bList and one subset of sets
        // 3.a. generate Equation that contains matrix a and vector b
        // 3.b. Solve with Gaussian Elimination
        //              3.b.i. 
        // 3.c. Check constraint validation
        // ==================================================================
        for(int index =0; index<sets.size(); index++){
        // -------------- Step 3.a. rebuild new Equations -------------------
        // ------------------------------------------------------------------   

        List<Integer> subset = sets.get(index);
        Equation myEQ = buildEquations(AList, bList, subset);
        
        // -------------- Step 3.b. Solve with Gaussian Elimination ---------
        // ------------------------------------------------------------------ 
        //System.out.println("0 . Matrix A: " + Arrays.deepToString(myEQ.a.toArray()) + " b: " + Arrays.toString(myEQ.b.toArray()));

        double[] soln = gaussianElimination(myEQ);
        SolnList.add(soln);
        
        
        
        // ------- Step 3.c. Check the Solution for Constraint Violation ----
        // ------------------------------------------------------------------ 
        boolean goodSoln = checkSolution(soln, AList, bList);
        Valid.add(goodSoln);
        
        
        
        double OBJ = 0.0;
        if(goodSoln)for(int i=0; i<c.length; i++) {
            if(Double.isInfinite(soln[i])&&c[i]!=0.0) return 1;
            OBJ += soln[i]*c[i];
        }
        OBJList.add(OBJ);
        
        
        }
        
        // ============== Step 4. Find Maxiumu Solution =====================
        // ==================================================================
        
        
        for(int iter = 0; iter<OBJList.size(); iter++){
//            System.out.println("Soln: "+Arrays.toString(SolnList.get(iter)));
//            System.out.println("Good? "+Valid.get(iter));
//            System.out.println("OBJ: "+OBJList.get(iter));
            if(Valid.get(iter)){
                if(maxOBJ <= OBJList.get(iter)){
                    noSoln = false;
                    maxOBJ = OBJList.get(iter);
                    maxIndex = iter;
                    
                    
                }
            }
        }
        
        if(noSoln) return -1;
        
        
        System.arraycopy(SolnList.get(maxIndex), 0, x, 0, SolnList.get(maxIndex).length);
        
        
        
        return 0;
    }
      
    static void preProcess(double[][] A, double[] b, List<List> AList, List<Double> bList){
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
        

    }
    
    static boolean checkSolution(double[] solution, List<List> A, List<Double> b) {
//        
//        double[][] A = myNewEQ.a;
//        double[] b = myNewEQ.b;
//        
//        for (int i = 0; i < A.length; i++) {
//            double[] inequality = A[i];
//            if (Arrays.stream(multiplyArrays(solution, inequality)).sum() > b[i])
//                return false;
//        }
        for (int constraint = 0; constraint<A.size(); constraint++){
            double sum = 0.0;
            for(int variable = 0; variable<A.get(0).size(); variable++){
                if(Double.isInfinite(solution[variable])){
                    //System.out.println("a infinite x");
                    if((double) A.get(constraint).get(variable) > 0){
                        //System.out.println("a non-zero coefficient");
                        return false;
                    }
                    sum += 0.0;
                    
                } else {
                    sum += solution[variable]*(double) A.get(constraint).get(variable);
                }
                    
            }
            if(sum > b.get(constraint)) return false;
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
                subset.add(iter);
            }
            Sets.add(subset);
            return Sets;
        } else {

            for (int i = 0; i < (1 << n); i++) {
                ArrayList<Integer> subset = new ArrayList<>();
                for (int j = 0; j < n; j++) {
                    if (((i >> j) & 1) == 1) {
                        subset.add(j);
                    } 
                }

                if (subset.size() == m) Sets.add(subset);                              
            }
            return Sets;
        }
    }
    
    static Equation buildEquations(List<List> AList, List<Double> bList, List<Integer> subset){
        List<ArrayList> AE = new ArrayList(); 
        List<Double> bE = new ArrayList(); 
        
        for(Integer index: subset){
                AE.add(new ArrayList((AList.get(index))));
                bE.add(bList.get(index));
        }

       
        Equation myEq = new Equation(AE,bE);
        
        return myEq;
       
    }
    
    static double[] gaussianElimination(Equation myEQ){
        
        
        double[] x = new double[myEQ.b.size()];
        Arrays.fill(x, 0.0);
        // -----------------------------------
        // STEP A: Select Pivot Row and Move it to the Top of the Un-processed Rows
        // -----------------------------------    
        
        
        for(int col = 0; col<myEQ.a.size(); col++){
            boolean skip = false;
            Position Pivot = new Position(-1,-1);
            for(int row = col; row<myEQ.a.size(); row++){
                if((double) myEQ.a.get(row).get(col) != 0.0) {
                        Pivot.row = row;
                        break; //this row is selected (pivoted) 
                }
            }
            
            // Parallel Condition
            if(Pivot.row == -1) {
//                double[] x = new double[myEQ.b.size()];
//                Arrays.fill(x, Double.POSITIVE_INFINITY);
//                return x;
                skip = true;
                x[col] = Double.POSITIVE_INFINITY;
                
            }
            
            if(!skip){
            ArrayList<Double> processedRow = new ArrayList<>(myEQ.a.get(Pivot.row)); 
            double temp_b = myEQ.b.get(Pivot.row);
            myEQ.a.remove(Pivot.row);
            myEQ.a.add(col, processedRow);
            myEQ.b.remove(Pivot.row);
            myEQ.b.add(col, temp_b);
            
            
        // -----------------------------------
        // STEP B: normalize that pivot element to 1
        // -----------------------------------  
            Pivot = new Position(col,col,(double) myEQ.a.get(col).get(col));
            
            for(int c=0; c<myEQ.a.size(); c++){
                double tempElement = (double) myEQ.a.get(Pivot.row).get(c);
                myEQ.a.get(Pivot.row).set(c, tempElement/Pivot.value);
            }
                myEQ.b.set(Pivot.row, (double)myEQ.b.get(Pivot.row)/Pivot.value);
            //System.out.println("1 . Matrix A: " + Arrays.deepToString(myEQ.a.toArray()) + " b: " + Arrays.toString(myEQ.b.toArray()));

        // -----------------------------------
        // STEP C: Process Other Rows
        // -----------------------------------  
            for(int oRow = 0; oRow <myEQ.a.size(); oRow ++){if(oRow != Pivot.row){
                double ratio = (double) myEQ.a.get(oRow).get(Pivot.column);
                for(int c=0; c<myEQ.a.size(); c++){
                double tempElement = (double) myEQ.a.get(oRow).get(c);
                    myEQ.a.get(oRow).set(c, tempElement - ratio*(double)myEQ.a.get(Pivot.row).get(c) );
                }
                    myEQ.b.set(oRow, (double) myEQ.b.get(oRow) - ratio*(double)myEQ.b.get(Pivot.row));
            }}
            //System.out.println("2 . Matrix A: " + Arrays.deepToString(myEQ.a.toArray()) + " b: " + Arrays.toString(myEQ.b.toArray()));
        }
        }
        

        // -----------------------------------
        // STEP D: Solve by back substitution
        // -----------------------------------

            
            
            for(int row = myEQ.a.size()-1; row>=0; row-- ){
                double oSum = 0.0;
                for(int col=0; col<myEQ.b.size(); col++){
                    if(col != row){
                        if((double) myEQ.a.get(row).get(col)!= 0){
                            oSum += (double) myEQ.a.get(row).get(col) * x[col];
                        } else {
                            oSum += 0;
                        }
                    }
                }
                //System.out.println("3. Matrix A: " + Arrays.deepToString(myEQ.a.toArray()) + " b: " + Arrays.toString(myEQ.b.toArray()));
                //System.out.println("oSUM: "+oSum+" coefficient: "+(double)myEQ.a.get(row).get(row)+" row: "+ row +" b: "+ myEQ.b.get(row));
                x[row] = (myEQ.b.get(row)- oSum)/(double)myEQ.a.get(row).get(row);
                //System.out.println(x[row]);
            }

        
        return x;
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
