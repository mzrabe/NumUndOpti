package org.mzrabe.methaheuristic;

import java.util.List;

import org.mzrabe.lina.Function;

/**
 * @author ,  
 */
public class SimulatedAnnealing
{
	/**
	 */
	public double[] solution;
	/** */
	public double bestValue;
	/**
	 * the maximal global number of iterations
	 */
	int maxGlobalIter = 1000000;
	/**
	 * the maximal local number of iterations, this means after how many
	 * iterations the algorithm stops if a change of configuration brings no
	 * benefit.
	 */
	int maxLokalIter = 2000;
	/**
	 * the function which checks the quality or fitness of the configuration
	 * (this is the function which generate the dimension to minimize)
	 */
	public Function func;
	/**
	 * the function which allows the simulated annealing algorithm a worser
	 * system configuration as in the last iteration step (known as temperature
	 * function)
	 */
	public Function tempFunc;
	/**
	 * the start temperature value
	 */
	public Double tempValue = 2000.;
	/**
	 * the name of the algorithms
	 */
	public final String ALGO_NAME = "Simulated Annealing";
	
	/**
	 * 
	 */
	public SimulatedAnnealing(Function f)
	{
	}
	
	public double[] findMin(double[] x, double[] c)
	{
		try
		{
			bestValue = func.getValue(x, c);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		
		
		
		
		return null;
	}

}
