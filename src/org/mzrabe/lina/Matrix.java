package org.mzrabe.lina;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mzrabe.lina.Gauss.InfinitySolutionsException;
import org.mzrabe.lina.Gauss.NoSolutionException;


/**
 * 
 * @author Moritz Zahn, mzrabe@gmail.com
 */
public class Matrix {
	
	private static final Logger log = LogManager.getRootLogger();
	
	/**
	 * Get the rotation matrix x-axes of the y-z Area.
	 * @param alpha - the angel in radiant
	 * @return - the rotation matrix
	 */
	public static double[][] getRx(double alpha)
	{
		return new double[][]{{1,0,0},{0, Math.cos(alpha),-Math.sin(alpha)},{0,Math.sin(alpha), Math.cos(alpha)}};
	}
	
	/**
	 * Get the rotation matrix y-axes of the x-z Area.
	 * @param beta - the angel in radiant
	 * @return - the rotation matrix
	 */
	public static double[][] getRy(double beta)
	{
		return new double[][]{{Math.cos(beta),0.,Math.sin(beta)},{0.,1.,0.},{-Math.sin(beta),0.,Math.cos(beta)}};
	}
	
	/**
	 * Get the rotation matrix z-axes of the x-y Area.
	 * @param gamma - the angel in radiant
	 * @return - the rotation matrix
	 */
	public static double[][] getRz(double gamma)
	{
		return new double[][]{{Math.cos(gamma),-Math.sin(gamma),0.},{Math.sin(gamma),Math.cos(gamma),0.},{0.,0.,1.}};
	}
	
	/**
	 * Get the inverse of the matrix A.
	 * @param A - the matrix
	 * @return the inverse matrix
	 */
	public static double[][] inverse(double[][] A)
	{
		//TODO use the efficiencies method depending on the size
		return inverseAT(A);
	}
	
	/**
	 * Get the inverse matrix of the matrix A. /TODO description 
	 * @param A -the matrix
	 * @return the inverse matrix
	 */
	public static double[][] inverseAT(double[][] A)
	{
		if (A[0].length != A.length)
			throw new IllegalArgumentException("The matrix has to be symetric.");
		if (det(A) == 0)
			throw new IllegalArgumentException("The matrix has no inverse, because the determinate is zero.");

		double[][] back = new double[A.length][A.length];

		/* make a copy */
		for (int i = 0; i < A.length; i++)
		{
			back[i] = Arrays.copyOf(A[i], back.length);
		}

		/** pivot element */
		double piv;
		int q = 0, p;

		/** index set */
		ArrayList<Integer> set = new ArrayList<>();

		for (int n = 0; n < back.length; n++)
		{
			/* find a pivot element which is not zero */
			p = 0;
			while ((piv = back[p][q]) == 0.0 || set.contains(p))
			{
				p++;
			}
			/* save the indexes */
			set.add(p);
			// System.out.println("p=" + p + ", q= " + q);
			back[p][q] = 1 / piv;
			// /* switch pivot-row */
			for (int k = 0; k < A.length; k++)
			{
				if (k != q)
					back[p][k] = -back[p][k] / piv;
			}
			/* switch rest elements */
			for (int i = 0; i < A.length; i++)
			{
				for (int k = 0; k < A.length; k++)
				{
					if (i != p && k != q)
						back[i][k] = back[i][k] + back[i][q] * back[p][k];
				}
			}
			/* switch pivot column */
			for (int i = 0; i < A.length; i++)
			{
				if (i != p)
					back[i][q] = back[i][q] / piv;
			}
			q++;
		}
		// printMatix(back);
		/* resort the matrix */
		Integer[] toSort = set.toArray(new Integer[0]);

		for (int i = 0; i < toSort.length - 1; i++)
		{
			for (int j = i + 1; j < toSort.length; j++)
			{
				if (toSort[i] > toSort[j])
				{
					/* sort the index */
					int temp = toSort[i];
					toSort[i] = toSort[j];
					toSort[j] = temp;

					/* switch column i with column j */
					for (int row = 0; row < back.length; row++)
					{
						double d = back[row][i];
						back[row][i] = back[row][j];
						back[row][j] = d;
					}
					/* switch row n with row m */
					int n = toSort.length - 1 - i;
					int m = toSort.length - 1 - j;
					for (int column = 0; column < back.length; column++)
					{
						double d = back[n][column];
						back[n][column] = back[m][column];
						back[m][column] = d;
					}
				}
			}
		}
		// System.out.println(Arrays.toString(toSort));
		// System.out.println("inv A by AT");
		// printMatix(back);
		return back;
	}
	
	/**
	 * Get the unit matrix of the given dimension.
	 * @param n - the dimension of the matrix
	 * @return - the unit matrix
	 */
	public static double[][] getE(int n)
	{
		double[][] back = new double[n][n];
		for(int i=0;i<n;i++)
		{
			back[i][i]=1;
		}
		return back;
	}
	
	/**
	 * Get the inverse matrix based on the Gauss-Jordan-Algorithms.
	 * @param matrix - the origin matrix to inverse
	 * @return - the inverse matrix
	 * @deprecated - FIXME bad result
	 */
	public static double[][] inverseGJ(double[][] matrix)
	{
		if (matrix[0].length != matrix.length)
			throw new IllegalArgumentException("The matrix has to be symetric.");
		
		double[][] A = new double[matrix.length][matrix.length];
		for(int i=0;i<A.length;i++)
		{
			A[i] = Arrays.copyOf(matrix[i], matrix.length);
		}

		double[][] B = new double[A.length][A.length];

		for (int i = 0; i < A.length; i++)
		{
			B[i][i] = 1;
		}
		
		/* make all entries below the diagonal to zero */
		for (int i = 0; i < A.length; i++)
		{
			 System.out.println("-------");
			// printMatix(A);
			/* find max absolute value in this column (pivot element) */
			int maxIDX = i;
			for (int j = 1; j < A.length; j++)
			{
				if (Math.abs(A[j][i]) > Math.abs(A[maxIDX][i]))
					maxIDX = j;
			}

			if (maxIDX != i)
			{
				/* switch row i with row maxIDX */
				double Atemp,Btemp;
				for (int j = 0; j < A.length; j++)
				{
					Atemp = A[i][j];
					A[i][j] = A[maxIDX][j];
					A[maxIDX][j] = Atemp;
					
					Btemp = B[i][j];
					B[i][j] = B[maxIDX][j];
					B[maxIDX][j] = Btemp;
				}
			}
			System.out.println("A");
			 printMatix(A);

			// TODO catch zero

			/* make all entries below aii to zero */

			for (int j = i + 1; j < A.length; j++)
			{
				/* calculate the quotient lij */
				A[j][i] = A[j][i] / A[i][i];
				/* calculate the reduced equations system */
				for (int l = 0; l < A.length; l++)
				{
					if(l>=i+1)
						A[j][l] = A[j][l] - A[j][i] * A[i][l];
					
					B[j][l] = B[j][l] - A[j][i] * B[i][l];
				}
				System.out.println("A");
				 printMatix(A);
			}
		}
		System.out.println("B");
		printMatix(B);
		
		/* make all entries over the diagonal to zero */
		System.out.println("******");
		for(int i=A.length-1;i>=0;i--)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				/* calculate the quotient uij */
				A[j][i] = A[j][i] / A[i][i];
				
				for (int l = A.length-1; l >= 0; l--)
				{
					if(l<=i-1 && l>j)
						A[j][l] = A[j][l] - A[j][i] * A[i][l];
					
					B[j][l] = B[j][l] - A[j][i] * B[i][l];
				}
			}
			System.out.println("A");
			printMatix(A);
			System.out.println("B");
			printMatix(B);
		}
		
		
		System.out.println("-------\nB");
		printMatix(B);
		
		for(int i=0;i<A.length;i++)
		{
			for(int j=0;j<A.length;j++)
			{
				B[i][j] /= A[i][i];
			}
			A[i][i]/=A[i][i];
		}
		
		System.out.println("A");
		printMatix(A);
		System.out.println("B");
		printMatix(B);
		
		System.out.println("A*");
		printMatix( multi(matrix, B));
		return B;
		
	}
	
	/**
	 * Calculate the determinate of the matrix A. The matrix has to be quadratic.
	 * @param 	A - the matrix
	 * @return the value of the det(A)
	 */
	public static double det(double[][] A)
	{
		
		if(A[0].length != A.length)
			throw new IllegalArgumentException("The matrix has to be quadratic.");
		
		double[][] copy = new double[A.length][A.length];
		for(int i=0;i<copy.length;i++)
		{
			copy[i]=Arrays.copyOf(A[i], copy.length);
		}
		
		double back = 1;
		
//		if(A.length == 2)
//			back+= A[0][0] * A[1][1] - A[0][1] * A[1][0];
//		else
//		{
			for(int i=0;i<copy.length-2;i++)
			{
//				System.out.println("-------");
//				printMatix(A);
				/* find max absolute value in this column (pivot element) */
				int maxIDX = i;
				for(int j=1;j<copy.length;j++)
				{
					if(Math.abs(copy[j][i]) > Math.abs(copy[maxIDX][i]))
						maxIDX = j;
				}
				
				if(maxIDX != i)
				{
					/* switch row i with row maxIDX */
					double temp;
					for(int j=0; j<copy.length;j++)
					{
						temp = copy[i][j];
						copy[i][j] = copy[maxIDX][j];
						copy[maxIDX][j] = temp;
					}
					back*=-1;
				}
				
//				printMatix(A);
				
				double aii = copy[i][i];
				//TODO catch zero
				back*=aii;
				
				/* make all entries below aii to zero */
				
				for(int j = i+1;j<copy.length;j++)
				{
					/* calculate the quotient lij */
					copy[j][i] = copy[j][i]/aii;
					/* calculate the reduced equations system */
					for(int l=i+1;l<copy.length;l++)
					{
						copy[j][l] = copy[j][l] - copy[j][i]*copy[i][l];
					}
				}
//				printMatix(A);
			}
			
			int n = copy.length-1;
			
			back*= copy[n-1][n-1] * copy[n][n] - copy[n-1][n] * copy[n][n-1];
//		}
//			System.out.println("det(A)= " + back);
		
		return back;
	}
	
	
	
	/**
	 * @param A - the left matrix[Rows][Columns]
	 * @param B - the right matrix[Rows][Columns]
	 * @return - the multiplication of the matrixes[Rows of A][Columns of B]
	 * @throws IllegalArgumentException - if the number of columns of A is not equal number of rows of B
	 */
	public static double[][] multi(double[][] A, double[][] B)
	{
			if(A[0].length != B.length)
				throw new IllegalArgumentException("Could not multiply matrix A["+A.length+"]["+A[0].length+"] * B["+B.length+"]["+B[0].length+"] because of the dimension.");
			
			double[][] C = new double[A.length][B[0].length];
			
			
			for(int i = 0;i<C.length;i++)
			{
				for(int j = 0;j<C[i].length;j++)
				{
					for(int k = 0;k<B.length;k++)
					{
						C[i][j]+= A[i][k] * B[k][j];
					}
				}
			}
			
			return C;
	}
	
	/**
	 * @param A - the left matrix[Rows][Columns]
	 * @param v - the right vector[Rows]
	 * @return - the multiplication
	 * @throws IllegalArgumentException - if the number of columns of A is not equal number of rows of v
	 */
	public static double[] multi(double[][] A, double[] v)
	{
		if(A.length == 0)
			return new double[0];
		
			if(A[0].length != v.length)
				throw new IllegalArgumentException("Could not multiply matrix A["+A.length+"]["+A[0].length+"] * v["+v.length+"] because of the dimension.");
			
			double[] C = new double[A.length];
			
			
			for(int i = 0;i<C.length;i++)
			{
				for(int j = 0;j<A[i].length;j++)
				{
					C[i]+= A[i][j] * v[j];
				}
			}
			
			return C;
	}
	
	/**
	 * @param A - the left matrix
	 * @param B - the right matrix
	 * @return - the addition of the matrix
	 * @throws IllegalArgumentException - if the number of columns of A is not equal number of rows of B
	 */
	public static double[][] addition(double[][] A, double[][] B)
	{	
		if(A.length != B.length && A[0].length != B[0].length)
			throw new IllegalArgumentException("Could not addition matrix A["+A.length+"]["+A[0].length+"] * B["+B.length+"]["+B[0].length+"] because of the dimension.");

		for(int i = 0; i<A.length;i++)
		{
			for(int j = 0; j<A[i].length;j++){
				A[i][j]+=B[i][j];
			}
		}
		return A;
	}
	
	
	/**
	 * @param A - the left matrix
	 * @param B - the right matrix
	 * @return - the subtraction of the matrix
	 * @throws IllegalArgumentException - if the number of columns of A is not equal number of rows of B
	 */
	public static double[][] subtract(double[][] A, double[][] B)
	{	
		if(A.length != B.length && A[0].length != B[0].length)
			throw new IllegalArgumentException("Could not addition matrix A["+A.length+"]["+A[0].length+"] * B["+B.length+"]["+B[0].length+"] because of the dimension.");

		for(int i = 0; i<A.length;i++)
		{
			for(int j = 0; j<A[i].length;j++){
				A[i][j]-=B[i][j];
			}
		}
		return A;
	}
	
	/**
	 * @param A - the matrix to transpose
	 * @return the transposed matrix A
	 */
	public static double[][] trans(double[][] A)
	{
		double[][] trans = new double[A[0].length][A.length];
		
		for(int i=0;i<A.length;i++)
		{
			for(int j=0;j<A[i].length;j++)
			{
				trans[j][i] = A[i][j];
			}
		}
		
		return trans;
	}
	
	/**
	 * @param A - matrix to print
	 */
	public static void printMatix(double[][] A)
	{	

		StringBuilder sb = new StringBuilder();
		for (int i=0;i<A.length;i++)
		{
			sb.append("|");
			for (int j=0;j<A[0].length-1;j++)
			{
				sb.append(String.format(Locale.ENGLISH,"%9f, ", A[i][j]));
			}
			sb.append(String.format(Locale.ENGLISH,"%9f |\n", A[i][A[0].length-1]));
		}
		System.out.println(sb.toString());
	}
	
	/**
	 * @param A - matrix to print
	 * @return - matrix as String
	 */
	public static String matixAsString(double[][] A)
	{	

		StringBuilder sb = new StringBuilder();
		for (int i=0;i<A.length;i++)
		{
			sb.append("|");
			for (int j=0;j<A[0].length-1;j++)
			{
				sb.append(String.format(Locale.ENGLISH,"%9f, ", A[i][j]));
			}
			sb.append(String.format(Locale.ENGLISH,"%9f |\n", A[i][A[0].length-1]));
		}
		return sb.toString();
	}
	
	/**
	 * Solve the linear equation system with the Gauss elimination and the relative column maximum strategy.
	 * Source: Schwarz2006 - Numerische Mathematik
	 * @param mat - the matrix
	 * @param b - the solution vector
	 * @return - the solution vector x
	 */
	public static double[] solveRCMS(double[][] mat, double[] b)
	{
		
		double[][] A = new double[mat.length][];
		for(int i=0;i<mat.length;i++)
		{
			A[i] = Arrays.copyOf(mat[i],mat[i].length);
		}
		
		log.entry();
		double max = 0;
		int[] p = new int[A.length];
		double temp;
		
		
		
		/* Gauss elimination */
		for(int k=0; k<A.length-1; k++)
		{
			max = 0;
			p[k] = 0;
			
			for(int i=k;i<A.length;i++)
			{
				double s = 0;
				for(int j = k;j<A.length;j++)
				{
					s = s + Math.abs(A[i][j]);
				}
				double q = Math.abs(A[i][k])/s;
				if(q > max)
				{
					max = q;
					p[k] = i;
				}
			}
			
			if(max == 0)
			{
				System.out.println("Cannot solve linear equation system, row is zero");
				return null;
			}
			if(p[k] != k)
			{
				
				for(int j = 0; j<A.length;j++)
				{
					temp = A[k][j];
					A[k][j] = A[p[k]][j];
					A[p[k]][j] = temp;
				}
			}
			for( int i = k+1;i<A.length;i++)
			{
				A[i][k] = A[i][k]/A[k][k];
				for(int j = k+1; j<A.length;j++)
				{
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
		
		/* switch the vector b */
		for(int k = 0; k<A.length-1; k++)
		{
			if(p[k] != k)
			{
				temp =b[k];
				b[k] = b[p[k]];
				b[p[k]] = temp;
			}
		}
		
		/* forward substitution */
		for( int i = 0; i<A.length;i++)
		{
			for(int j=0;j<i-1;j++)
			{
				b[i] -= A[i][j] * b[j];
			}
		}
		
		/* backward substitution */
		for( int i = A.length-1;i>=0;i--)
		{
			temp = b[i];
			for(int k = i+1;k<A.length;k++)
			{
				temp -= A[i][k] * b[k];
			}
			b[i] = temp/A[i][i];
		}
		log.exit();
		return b;
	}
	
	/**
	 * Solve the linear equation system with the Gauss elimination and the row maximum strategy.
	 * Source: Schwarz2006 - Numerische Mathematik
	 * @param mat - the matrix
	 * @param b - the solution vector
	 * @return - the solution vector x
	 */
	public static double[] solveRMS(double[][] mat, double[] b)
	{
		/* LR-Zerlegung */
		
		double[][] A = new double[mat.length][];
		for(int i=0;i<mat.length;i++)
		{
			A[i] = Arrays.copyOf(mat[i],mat[i].length);
		}
		
		int[] p = new int[A.length-1];
		
		for(int i = 0; i<A.length-1;i++)
		{	
			
			double max = Math.abs(A[i][i]);
			int maxID = i;
			
			for(int j = i+1; j<A.length;j++)
			{
				if(max < Math.abs(A[j][i]))
				{
					maxID = j;
					max = Math.abs(A[j][i]);
				}
			}
			
			p[i] = maxID;
			
			if(p[i] != i)
			{
				double temp;
				for(int j=0;j<A.length;j++)
				{
					temp = A[p[i]][j];
					A[p[i]][j] = A[i][j];
					A[i][j] = temp;
				}
				
				temp = b[p[i]];
				b[p[i]] = b[i];
				b[i] = temp;
			}
			
			
			// spalte
			for(int k = i+1; k<A.length;k++)
			{
				// zeile
				A[k][i] = A[k][i] / A[i][i];
				
				for(int j=i+1; j<A[k].length;j++)
				{
					// spalte
					A[k][j] =  A[k][j] - A[i][j] * A[k][i];
					
				}
				b[k] =  b[k] - b[i] * A[k][i];
			}
		}
		
		for(int i=A.length-1;i>=0;i--)
		{
			double sum = 0;
			
			if(i != A.length)
			{
				for(int k=i+1;k<A.length;k++)
				{
					sum+= A[i][k]*b[k];
				}
			}
			b[i] = (b[i] - sum)/A[i][i];
		}
		
//		Vector.print(switchLog);
		return b;
	}
	
	/**
	 * Solve the linear equation system with the Chorlesky decomposition. Only possible if the matrix A is a symmetric positive definite matrix.
	 * Source: Schwarz2006 - Numerische Mathematik
	 * @param mat - the matrix
	 * @param b - the vector
	 * @return - the solution vector x
	 */
	public static double[] solveChorlesky(double[][] mat, double[] b)
	{
		log.entry();
		
		double[][] A = new double[mat.length][];
		for(int i=0;i<mat.length;i++)
		{
			A[i] = Arrays.copyOf(mat[i],mat[i].length);
		}
		
		for(int k=0;k<A.length;k++)
		{
			if(A[k][k] <= 0)
				log.info("Matrix is not positiv definit");

			A[k][k] = Math.sqrt(A[k][k]);
			
			for(int i = k+1; i<A.length; i++)
			{
				A[i][k] = A[i][k]/A[k][k];
				for(int j=k+1;j<=i;j++)
				{
					A[i][j] = A[i][j] - A[i][k] * A[j][k];
				}
			}
		}
		
		for(int i = 0;i<A.length;i++)
		{
			double s = b[i];
			for(int j = 0; j<=i-1;j++)
			{
				s = s - A[i][j] * b[j];
			}
			b[i] = s/A[i][i];
		}
		
		for(int i = A.length-1; i>=0; i--)
		{
			double s = b[i];
			for(int k = i+1;k<A.length;k++)
			{
				s = s - A[k][i] * b[k];
			}
			b[i] = s/A[i][i];
		}
		
		log.exit();
		return b;
	}
	
	/**
	 * Function for finding rank of coefficient vector matrix. This method use the library of JAMA.
	 * @param mat - matrix
	 * @param vector - vector
	 * @return - rank of the given coefficient vector matrix
	 */
	public static int rankOfMatrix(double mat[][], double[] vector)
	{
		if(mat.length != vector.length)
			throw new IllegalArgumentException("Matrix and Vector must have the same number of lines.");
		
		double[][] A = new double[mat.length][mat[0].length+1];
		
		for(int i = 0;i<A.length;i++)
		{
			for(int j = 0; j<A[0].length;j++)
			{
				if(j==A[0].length-1)
					A[i][j] = vector[i];
				else
					A[i][j] = mat[i][j];
			}
		}
		
		return rankOfMatrix(A);
	}

	/**
	 * Function for finding rank of matrix. This method use the library of JAMA.
	 * @param mat - matrix
	 * @return - rank of the given matrix
	 */
	public static int rankOfMatrix(double mat[][])
	{
		return new Jama.Matrix(mat).rank();
	}
	
	/**
	 * Solve the linear system of equations taking into account the known variables.
	 * @param A - coefficient maxtrix
	 * @param b - right hand vector of the linear system
	 * @param x - solution vector, unknown variables setted as Double.NaN otherwise the known value. <br><b>NOTE! The solution will save into this double array.</b>
	 * @return - the solution of the linear system of equations
	 * @throws InfinitySolutionsException 
	 * @throws NoSolutionException 
	 */
	public static double[] solveLinearSystem(double[][] A, double[] b, double[] x) throws NoSolutionException, InfinitySolutionsException
	{	
		/* at first count the known variables of x */
		
		int known_x = 0;
		
		for(int i = 0; i<x.length;i++)
		{
			if(!Double.isNaN(x[i]))
				known_x++;
		}
		
		/* check if it is possible to solve the linear system of equation */
		
		if(A.length < x.length - known_x)
			throw new IllegalArgumentException("The number of equations ("+A.length+") is not greather than the number of unknown variables ("+(x.length - known_x)+"). Can't solve the linear system of equation.");
			
		/* Now reduce the linear equation system according to the
		 * fix property nodes.
		 */
		
		int i_ = 0,j_ = 0;
		double[][] A_ = new double[x.length][x.length - known_x];
		double[] b_ = new double[x.length];
		
		for(int i = 0;i<A.length;i++)
		{
//			if(i==A_.length)
//				break;
//			if(!Double.isNaN(x[i]))
//				continue;
			b_[i_] = b[i];
			j_ = 0;
			for(int j = 0;j<A[0].length;j++)
			{
				if(!Double.isNaN(x[j]))
				{
					b_[i_] -= A[i][j] * x[j];
					continue;
				}
				
				A_[i_][j_] = A[i][j];
				j_++;
			}
			i_++;
		}
		
		/* count zero lines */
		ArrayList<Integer> zero = new ArrayList<>();
		for(int i = 0; i<A_.length;i++)
		{
			if(b_[i] == 0)
			{
				double sum = 0;
				for(int j = 0;j<A_[0].length;j++)
				{
					sum+=Math.abs(A_[i][j]);
				}
				if(sum == 0.0)
					zero.add(i);
			}
		}
		
		int rankAb = rankOfMatrix(A_, b_);
		int rankA = rankOfMatrix(A_); 
		
		if(rankAb != rankA)
		{
//			Matrix.printMatix(A_);
//			Vector.print(b_);
			//throw new Gauss().new NoSolutionException("No solution. rank(Ab) = " + rankAb + ", rank(A) = " + rankA);
		}
//		double[][] A_final = new double[rankA][A_[0].length];
//		double[]   b_final = new double[rankA];
//		
//		i_ = 0;
//		for(int i = 0; i<A_.length;i++)
//		{
//			if(i_ == A_final.length)
//				break;
//			
//			if(!zero.contains(i))
//			{
//				for(int j = 0;j<A_[0].length;j++)
//				{
//					A_final[i_][j]=A_[i][j];
//				}
//				
//				b_final[i_] = b_[i];
//				i_++;
//			}
//		}
		
		double[] x_ = new Jama.Matrix(A_).solve(new Jama.Matrix(b_,b_.length)).getColumnPackedCopy();//Gauss.getSolution(A_final, b_final, false);
		
		i_ = 0;
		for(int i = 0;i<x.length;i++)
		{
			if(Double.isNaN(x[i]))
			{
				x[i] = x_[i_];
				i_++;
			}
		}
		
		return x;
	}
}
