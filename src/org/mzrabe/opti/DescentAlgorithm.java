package org.mzrabe.opti;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Vector;

/**
 * 
 * @author Moritz Zahn <mzrabe@gmail.com>
 */
public class DescentAlgorithm {
	

	private static final Logger log = LogManager.getRootLogger();
	
	private static double toleranc = 1e-5;
	private static int maxIterations = (int) 1e3;
	private static StepStrategy step = new ArmijoStep();
	public static ArrayList<double[]> solutionHistory = new ArrayList<>();
	
	public static double[] localNewton(Function f, double[] x) throws Exception
	{
		log.info("Start local newton optimization ...");
		double[] solution = Arrays.copyOf(x,x.length);
		double[] d = new double[x.length];
		int numIterations = 1;
		solutionHistory.clear();
		solutionHistory.add(solution);
		
		
		
		
		do
		{
			d = Gauss.getSolution(f.getHf(solution), f.grad(solution,true), false);
			solution = Vector.sum(solution, d);
			numIterations++;
			solutionHistory.add(solution);
		}
		while(Vector.twoNorm(d) > toleranc && numIterations < maxIterations);
		
		log.info("Local newton algorithms found minimum after "+numIterations+" iterations at " + Arrays.toString(solution) + "\n");
		
		return solution;
	}
	
	public static double[] simpleLocalNewton(Function f, double[] x) throws Exception
	{
		log.info("Start simple newton optimization ...");
		double[] solution = Arrays.copyOf(x,x.length);
		double[] d = new double[x.length];
		
		int numIterations = 1;
		
		
		do
		{
			d = Gauss.getSolution(f.getHf(x), f.grad(solution,true), false);
			solution = Vector.sum(solution, d);
			numIterations++;
		}
		while(Vector.twoNorm(d) > toleranc && numIterations < maxIterations);
		
		log.info("Simple local newton algorithms found minimum after "+numIterations+" iterations at " + Arrays.toString(solution)+ "\n");
		
		return solution;
	}
	
	public static double[] steamedLocalNewton(Function f, double[] x) throws Exception
	{
		log.info("Start steamed newton optimization ...");
		double[] solution = Arrays.copyOf(x,x.length);
		double[] d = new double[x.length];
		int numIterations = 1;
		solutionHistory.clear();
		solutionHistory.add(solution);
		
		
		do
		{
			d = Gauss.getSolution(f.getHf(solution), f.grad(solution,true), false);
			log.debug("d = " +  Arrays.toString(d));
			d = Vector.multiScalar(d, step.getAlpha(f,d,solution));
			log.debug("d * alpha= " +  Arrays.toString(d));
			solution = Vector.sum(solution, d);
			log.debug("x + d * alpha= " +  Arrays.toString(solution));
			numIterations++;
			solutionHistory.add(solution);
		}
		while(Vector.twoNorm(d) > toleranc && numIterations < maxIterations);
		
		log.info("Steamed local newton algorithms found minimum after "+numIterations+" iterations at " + Arrays.toString(solution)+ "\n");
		
		return solution;
	}
	
	public static double[] gradientDescent(Function f, double[] x, double[] c ) throws Exception
	{
		log.info("Start gradient descent optimization ...");
		double[] solution = Arrays.copyOf(x, x.length);
		double[] d = f.grad(solution, true);
		double alpha = 1;
		int numIterations = 0;
		solutionHistory.clear();
		solutionHistory.add(solution);
		
		while(Vector.twoNorm(d) > toleranc && numIterations < maxIterations)
		{	
			numIterations++;
			d = f.grad(solution, true);
			alpha = step.getAlpha(f, d, solution, c);
			d = Vector.multiScalar(d, alpha);
			solution = Vector.sum(solution, d);
			solutionHistory.add(solution);
		}
		
		log.info("Gradient descent algorithms found minimum after "+numIterations+" iterations at " + Arrays.toString(solution)+ "\n");
		
		return solution;
		
		
	}

	public static double[] globalNewton(Function f, double[] x, double[] c, double delta, double roh, boolean steamedLocalNewton) throws Exception
	{
		log.info("Start global newton optimization ...");
		double[] solution = Arrays.copyOf(x, x.length);
		double[] d = new double[solution.length];
		double alpha = 1;
		int numIterations = 0;
		solutionHistory.clear();
		solutionHistory.add(solution);
		
		do
		{	
			numIterations++;
			
			d = Gauss.getSolution(f.getHf(x), f.grad(solution,true), false);
			
			if(d == null || Vector.scalarProduct(f.grad(solution, false),d) > Math.pow(-delta*(Vector.twoNorm(d)), roh))
			{
				d = f.grad(solution, true);
				alpha = step.getAlpha(f, d, solution, c);
				d = Vector.multiScalar(d, alpha);
			}
			else
			{
				if(steamedLocalNewton)
				{
					alpha = step.getAlpha(f, d, solution, c);
					d = Vector.multiScalar(d, alpha);
				}
				
			}
			
			solution = Vector.sum(solution, d);
			solutionHistory.add(solution);
		}
		while(Vector.twoNorm(d) > toleranc && numIterations < maxIterations);
		
		log.info("Global newton algorithms found minimum after "+numIterations+" iterations at " + Arrays.toString(solution)+ "\n");
		
		return solution;
	}

}
