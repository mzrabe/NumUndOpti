package org.mzrabe.opti;

import java.util.Arrays;
import java.util.TreeMap;
import org.mzrabe.lina.Function;

import static org.mzrabe.lina.Vector.*;

/**
 * 
 * @author Moritz Zahn <mzrabe@gmail.com>
 *
 */
public class DownhillSimplex extends OptiAlgorithm{
	
	
	/**
	 * the first index of the points
	 */
	private final int FIRST_P = 0;
	/**
	 * the last index of the points
	 */
	private int LAST_P;
	/**
	 * the number of simplex points
	 */
	private int NPoints;
	/**
	 * the size of dimension
	 */
	private int dimesion;
	/**
	 * the points of the simplex
	 */
	private double[][] points;
	/**
	 * coefficients for the algorithm (reflection, expansion, contraction)
	 */
	double alpha, beta, gamma;
	double[] psi_, psiR;
	
	/**
	 * Initialize the downhill simplex algorithm of Nelder and Mead.
	 * 
	 * @param func
	 *            - the function where a minimum is searched
	 * @param tolerance
	 *            - the tolerance the function value
	 * @param maxNumIteration
	 *            - the maximal number of iterations
	 */
	public DownhillSimplex(Function func,double tolerance, int maxNumIteration)
	{
		this(func, tolerance, maxNumIteration, 1.5, 0.9);
	}
	
	
	/**
	 * Initialize the downhill simplex algorithm of Nelder and Mead.
	 * 
	 * @param func
	 *            - the function where a minimum is searched
	 * @param tolerance
	 *            - the tolerance the function value
	 * @param maxNumIteration
	 *            - the maximal number of iterations
	 * @param alpha
	 * 			  - coefficient for the reflection
	 * @param gamma
	 * 			  - coefficient for the contraction
	 */
	public DownhillSimplex(Function func,double tolerance, int maxNumIteration, double alpha, double gamma)
	{
		super(func, tolerance, maxNumIteration);

		
		this.alpha = alpha;
		this.beta = Math.max(alpha, 1);
		this.gamma = gamma;
		
		log.debug("number of points: " + NPoints);
		log.debug("dimension: " + dimesion);
		log.debug("tolerance: " + tolerance);
		log.debug("maxNumIterations: " + maxNumIterations);
		log.debug("alpha: " + alpha);
		log.debug("beta: " + beta);
		log.debug("gamma: " + gamma);
		
	}
	
	/**
	 * Contract the Simplex.
	 */
	private void shrinkSimplex()
	{
		log.debug("Shrink simplex ...");
		for(int i=1;i<NPoints;i++)
		{
			points[i] = multiScalar(sum(points[FIRST_P], points[i]), 0.5);
		}
	}
	
	public double[] getPsi_()
	{
		double[] sum = new double[dimesion];
		
		for(int i=0;i<NPoints-1;i++){
			sum = sum(sum, points[i]);
		}
		
		return multiScalar(sum, 1./dimesion);
	}
	
	@Override
	protected boolean isFinish()
	{
		return (calcAverage() < tolerance)  ? true : false;
	}
	
	private double calcAverage()
	{
		double[] funcValues = new double[NPoints];
		double sumFuncValue = 0;
		
		for(int i=0;i<NPoints;i++)
		{
			funcValues[i] = func.getValue(points[i]);
			sumFuncValue += funcValues[i];
		}
		sumFuncValue =sumFuncValue/NPoints;
		
		double sumQFuncValue = 0;
		for(int i=0;i<NPoints;i++)
		{
			sumQFuncValue += Math.pow(funcValues[i]-sumFuncValue, 2);
		}
		
		sumQFuncValue = sumQFuncValue/NPoints;
		
		log.debug("sumQFuncValue : "+sumQFuncValue);
		return sumQFuncValue;
	}
	
	/**
	 * Resort the points. They will sorted ascending depending on the function values at this point.
	 */
	public void sortPoints()
	{
		double[] f = new double[points.length];
		
		for(int i=0;i<points.length;i++)
		{
			f[i] = func.getValue(points[i]);
		}
		
		
		for(int i = 0;i<points.length-1;i++)
		{
			double min = f[i];
			int idx = i;
			
			for(int j=i+1;j<points.length;j++)
			{
				if(f[j]<min)
				{
					min = f[j];
					idx = j;
				}
			}
			if(idx != i)
			{
				double[] minPoint = points[idx];
				points[idx] = points[i];
				f[idx] = f[i];
				points[i] = minPoint;
				f[i] = min;
			}
			
			
			
			log.trace(f[idx] + " : " + Arrays.toString(points[i]));
		}
		
		log.trace(f[points.length-1] + " : " + Arrays.toString(points[points.length-1]));
	}
	
	/**
	 * @param x
	 *            - the start points for the downhill simplex algorithm, the
	 *            number of points must be greater by 1 then the dimension of
	 *            the points.
	 */
	@Override
	protected void initAlgorithms(double[] ... x) {
		this.NPoints = x.length;
		this.dimesion = x[FIRST_P].length;
		this.LAST_P = NPoints -1;
		
		/* the number of points must be greater than the dimension of the points by 1 */
		if(NPoints != (dimesion+1))
		{
			throw new NumberFormatException("The number of points must be greater than the dimension of the points by 1.");
		}
		
		points = new double[x.length][];
		for(int i=0;i<x.length;i++)
		{
			points[i] = x[i].clone();
		}
		
		sortPoints();
//		log.debug("start simplex : " + vectorsToString(points));
		solution = this.getPsi_()/*points[FIRST_P]*/;
	}


	@Override
	protected void algorithms() {
		
		
		psi_ = getPsi_();
		psiR = sum(psi_, multiScalar(minus(psi_, points[LAST_P]),alpha));
	
		if(func.getValue(points[FIRST_P]) <= func.getValue(psiR) && func.getValue(psiR) <= func.getValue(points[LAST_P-1]))
		{
			log.debug("reflection: points[FIRST_P] <= psiR <= points[LAST_P]");
			points[LAST_P] = psiR;
		}
		else if(func.getValue(psiR) < func.getValue(points[FIRST_P]))
		{
			log.debug("expansion: psiR < points[FIRST_P]");
			double[] psiE = sum(psi_, multiScalar(minus(psiR, psi_),beta));
			
			points[LAST_P] = func.getValue(psiE) < func.getValue(psiR) ?  psiE : psiR;
		}
		else if(func.getValue(psiR) > func.getValue(points[FIRST_P]))
		{
			log.debug("contraction: psiR > points[FIRST_P]");
			double[] psiC;
			
			if(func.getValue(psiR) >= func.getValue(points[LAST_P]))
			{
				log.debug("psiR >= points[LAST_P]");
				psiC = sum(psi_, multiScalar(minus(points[LAST_P], psi_),gamma));
			}
			else
			{
				log.debug("psiR < points[LAST_P]");
				psiC = sum(psi_, multiScalar(minus(psiR, psi_),gamma));
			}
			
			if(func.getValue(psiC) < Math.min(func.getValue(points[LAST_P]), func.getValue(psiR)))
			{
				points[LAST_P] = psiC;
			}
			else
			{
				shrinkSimplex();
			}
		}
		
		sortPoints();

		
	}
	
	/**
	 * Get the simplex of the solution.
	 * @return
	 */
	public double[][] getSolutionSimplex()
	{
		double[][] copy = new double[points.length][];
		for(int i=0;i<points.length;i++)
		{
			copy[i] = Arrays.copyOf(points[i],points[i].length);
		}
		return copy;
	}


	@Override
	public String getName() {
		return "Downhill Simplex";
	}


	@Override
	public double[] getSolution()
	{
		return getPsi_();
	}
	

}
