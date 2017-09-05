package org.mzrabe.zeropoint;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mzrabe.diffquotient.JacobianMatrix;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Matrix;
import org.mzrabe.lina.Vector;

/**
 * Implementation of the Newton algorithms
 * @author Moritz Zahn
 */
public class Newton
{
	/** precision of the solution */
	private static double precision = 1e-5;
	/** maximal number of iterations */
	private static int maxIter = 1000;
	/** iteration steps for logging */
	private static int logStep = 50;
//	/** the implemetation to get the Jacobian matrix */
//	private JacobianMatrix j = new JacobianMatrix();
	private static final Logger log = LogManager.getRootLogger();
	
	private static BufferedWriter bw;
	private static BufferedWriter bw1;
	
	
	/**
	 * @param f
	 * @param x0
	 * @param c
	 * @return
	 * @throws Exception 
	 */
	public static double[] newtonStep(Function[] f, double[] x0, double... c) throws Exception
	{
		double[][] j = JacobianMatrix.getJacobiMatrix(f, x0, c);
		double[] mf = getMinusF(f, x0, c);
//		double[] r = Gauss.getSolution(j, mf, false);
		double[] r = Matrix.solveRMS(j, mf);
		
		return r;
	}
	
	/**
	 * See {@link Newton#getSolution(Function[], double[], double...)}
	 * @param f - the function
	 * @param printSolution - flag if the solution should print into the log
	 * @param x - start point for the iteration
	 * @param c - constant values of f
	 * @return - the zero point of f
	 * @throws Exception 
	 */
	public static double getSolution(Function f,boolean printSolution, double x, double ... c) throws Exception
	{
		return getSolution(new Function[]{f}, printSolution, new double[]{x}, c)[0];
	}

	/**
	 * Approximate a solution of the equation system according to the newton
	 * method.<br>
	 * <br>
	 * 
	 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"> <msub>
	 * <mi>J</mi> <mi>f</mi> </msub> <mo stretchy="false">(</mo> <msup>
	 * <mi>x</mi> <mrow class="MJX-TeXAtom-ORD"> <mi>(n)</mi> </mrow> </msup>
	 * <mo stretchy="false">)</mo> <mo>&#x22C5;<!-- ⋅ --></mo>
	 * <mo stretchy="false">(</mo> <msup> <mi>x</mi>
	 * <mrow class="MJX-TeXAtom-ORD"> <mo stretchy="false">(</mo> <mi>n</mi>
	 * <mo>+</mo> <mn>1</mn> <mo stretchy="false">)</mo> </mrow> </msup>
	 * <mo>&#x2212;<!-- - --></mo> <msup> <mi>x</mi>
	 * <mrow class="MJX-TeXAtom-ORD"> <mo stretchy="false">(</mo> <mi>n</mi>
	 * <mo stretchy="false">)</mo> </mrow> </msup> <mo stretchy="false">)</mo>
	 * <mo>=</mo> <mi> &#x2212; f </mi> <mo>(</mo> <msup> <mi>x</mi>
	 * <mrow class="MJX-TeXAtom-ORD"> <mi>(n)</mi> </mrow> </msup>
	 * <mo stretchy="false">)</mo> </math>
	 * 
	 * Its manly used for a non linear equation system. This method use the
	 * Gauss-Jordan algorithms to solve the linear equation system of one
	 * iteration. If this linear equation system cannot solved null will
	 * returned. If the maximal number of iterations are reached also null is
	 * return, because you cannot say anything about the quality of the current
	 * approximation.
	 * 
	 * @param f
	 *            - a vectorial function or the equation system
	 * @param printSolution
	 *            - flag if the solution should print into the log
	 * @param x0
	 *            - the initial values for the Newton algorithms
	 * @param c
	 *            - constant values in the vectorial function (f)
	 * @return - a vector which solved the equation system
	 * @throws Exception 
	 */
	public static double[] getSolution(Function[] f, boolean printSolution, double[] x0, double ... c) throws Exception
	{
		int numIter = 1;
		double[][] j = JacobianMatrix.getJacobiMatrix(f, x0, c);
		double[] mf = getMinusF(f, x0, c);
//		double[] r = Gauss.getSolution(j, mf, false);
		double[] r = Matrix.solveRMS(j, mf);
		double[] rMin = Arrays.copyOf(r, r.length);
		double[] bestSolution = Arrays.copyOf(x0, x0.length);
		
//		try
//		{
//			bw = new BufferedWriter(new FileWriter(new File("log_func.txt")));
//			bw1 = new BufferedWriter(new FileWriter(new File("log_r.txt")));
//		}
//		catch (IOException e1)
//		{
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
		
//		writeLog(getMinusF(f, x0, c), r);
		
		if(r == null)
			return null;
		
		while(Vector.oneNorm(r) > precision)
		{
			if(numIter >= maxIter)
			{
				if(printSolution)
				{
					log.info("Maximal number of iteration are reached.");
					StringBuilder sb = new StringBuilder();
					sb.append("Found best solution at x = ");
					sb.append(Vector.asString(bestSolution));
					sb.append(" after ");
					sb.append(numIter);
					sb.append(" iterations.");
					
					log.info(sb.toString());
					
					sb.setLength(0);
					
					sb.append("r = ");
					sb.append(Vector.asString(rMin));
					
					log.info(sb.toString());
					
					sb.setLength(0);
					sb.append("f = ");
					
					double calcF[] = new double[f.length];
					
					for(int i=0;i<f.length-1;i++)
					{
						calcF[i] = f[i].getValue(x0, c);
					}
					
					sb.append(Vector.asString(calcF));
					
					log.info(sb.toString());
				}
				
				return bestSolution;
			}
			
			if((numIter%logStep) == 0 && printSolution)
			{
				log.info("----------------------------");
				StringBuilder sb = new StringBuilder();
				sb.append("x = ");
				sb.append(Vector.asString(x0));
				sb.append(" after ");
				sb.append(numIter);
				sb.append(" iterations.");
				
				log.info(sb.toString());
				
				sb.setLength(0);
				sb.append("r = ");
				sb.append(Vector.asString(r));
				
				log.info(sb.toString());
				
				
				sb.setLength(0);
				sb.append("f = ");
				
				double calcF[] = new double[f.length];
				
				for(int i=0;i<f.length-1;i++)
				{
					calcF[i] = f[i].getValue(x0, c);
				}
				
				sb.append(Vector.asString(calcF));
				
				log.info(sb.toString());
				
			}
//			
//			Vector.print(r);
//			Vector.print(x0);
//			Vector.print(getMinusF(f, x0, c));
//			System.out.println("--------------");
			
			Vector.add(x0, r);
			
			
			
			
			j = JacobianMatrix.getJacobiMatrix(f, x0, c);
			mf = getMinusF(f, x0, c);
//			r = Gauss.getSolution(j, mf, false);
			r = Matrix.solveRMS(j, mf);
			
			if( Vector.oneNorm(r) < Vector.oneNorm(rMin))
			{
				rMin = Arrays.copyOf(r, r.length);
				bestSolution = Arrays.copyOf(x0, x0.length);
			}
			
			
			
			if(r == null)
			{
				log.info("Linear equation system cannot solved. r is null");
				return null;
			}
			
//			writeLog(getMinusF(f, x0, c), r);
			
			numIter++;
		}
			
		if(printSolution)
		{
			StringBuilder sb = new StringBuilder();
			sb.append("Found solution at x = ");
			sb.append(Vector.asString(bestSolution));
			sb.append(" after ");
			sb.append(numIter);
			sb.append(" iterations.");
			
			log.debug(sb.toString());
			
			sb.setLength(0);
			
			sb.append("r = ");
			sb.append(Vector.asString(rMin));
			
			log.debug(sb.toString());
			
			sb.setLength(0);
			sb.append("f = ");
			
			double calcF[] = new double[f.length];
			
			for(int i=0;i<f.length-1;i++)
			{
				calcF[i] = f[i].getValue(x0, c);
			}
			
			sb.append(Vector.asString(calcF));
			
			log.debug(sb.toString());
		}
		
//		try
//		{
//			bw.flush();
//			bw.close();
//			bw1.flush();
//			bw1.close();
//		}
//		catch (IOException e)
//		{
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		
		return bestSolution;
	}
	
	public static double[] getMinusF(Function[] f, double[] x, double ... c) throws Exception
	{
		double[] mf = new double[f.length];
		
		for(int i=0;i<f.length;i++)
		{
//			try
//			{
			mf[i] = -f[i].getValue(x,c);
//			}
//			catch(IllegalArgumentException e)
//			{
//				System.out.println(String.format("Arguments of x=%s are not valide.", Vector.asString(x)));
//				return null;
//			}
		}
		
		return mf;
	}
	
	private static void writeLog(double[] mf, double[] r)
	{
		try
		{
			for(double d : mf)
			{
				bw.write(String.valueOf(-d)+",");
			}
			bw.write("\n");
			bw.flush();
			
			for(double d : r)
			{
				bw1.write(String.valueOf(Math.abs(d))+",");
			}
			bw1.write("\n");
			bw1.flush();
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Set the precision of the solution
	 * @param d - the precision, they has to be greater than 1.e-15
	 */
	public static void setPrecision(double d)
	{
		if(precision <= 1e-15)
			throw new IllegalArgumentException("The precision has to be greater than 1.0e-15.");
		
		Newton.precision = d;
	}
	
	/**
	 * @return - the currently set precision of the solution
	 */
	public static double getPrecision()
	{
		return Newton.precision;
	}
}
