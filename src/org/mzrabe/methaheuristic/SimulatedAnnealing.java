package org.mzrabe.methaheuristic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import org.mzrabe.exception.SetUpException;

public class SimulatedAnnealing<T, R extends Number & Comparable<R>> extends Thread
{
	/**
	 * the array with possible objects which order influence the fitness of the
	 * system
	 */
	List<T> objects;
	/**
	 * array of indexes which represent the configuration of the system
	 */
	public int[] optConf;
	/**
	 * the start configuration
	 */
	public int[] startConf;
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
	public Function<List<T>, R> fitFunc;
	/**
	 * the function which allows the simulated annealing algorithm a worser
	 * system configuration as in the last iteration step (known as temperature
	 * function)
	 */
	public Function<Integer, Double> tempFunc;
	/**
	 * the start temperature value
	 */
	public Double tempValue = 2000.;
	/**
	 * the data to plot
	 */
	public List<R> fitData;
	/**
	 * the name of the algorithms
	 */
	public final String ALGO_NAME = "Simulated Annealing";

	

	/**
	 * Constructor to initialize the simulated annealing algorithms.
	 * 
	 * @param objects
	 *            - the list of objects
	 * @param fitFunc
	 *            - the fit function makes a statement about goodness of the
	 *            configuration
	 */
	public SimulatedAnnealing(List<T> objects, Function<List<T>, R> fitFunc)
	{
		this.objects = objects;
		this.fitFunc = fitFunc;
		this.fitData = new ArrayList<>();
		this.setName(ALGO_NAME);
	}
	
	@Override
	public void run()
	{
		if (fitFunc == null)
		{
			System.out.println("No fitness function was set up. Cannot run optimization algorithm.");
			interrupt();
		}
		if (tempFunc == null)
			tempFunc = getDefaultTemperatureFunction();
		/* start the algorithms */
		
		/* generate a random start configuration */
		optConf = generateRandomPermutation(objects.size());
		startConf = Arrays.copyOf(optConf, optConf.length);

		/* inti first the yData with one value */
		this.fitData.add(getOptimizedValue());
		

		/* init the iterations counters */
		int numG = 0;
		int numL = 0;

		int[] si = Arrays.copyOf(optConf, optConf.length);
		int[] sj;
		// int[] tempConf;
		
		File f = new File("yData.txt");
		try
		{
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));

		
		do
		{
			/* save the best solution */
			if (fitFunc.apply(getConfiguration(si)).compareTo(getOptimizedValue()) < 0)
			{
				optConf = Arrays.copyOf(si, si.length);
			}

			/* change the configuration of si */
			sj = Arrays.copyOf(si, si.length);
			changeConfiguration(sj);

			R qi = fitFunc.apply(getConfiguration(si));
			R qj = fitFunc.apply(getConfiguration(sj));

			/* check if the new configuration is better as before */
			if (qj.compareTo(qi) < 0)
			{
				/* save the new configuration */
				si = Arrays.copyOf(sj, sj.length);
				numL = 0;
			}
			else if (Math.random() < Math.exp(-(qj.doubleValue() - qi.doubleValue()) / tempFunc.apply(numG)))
			{
				/* save the new configuration */
				si = Arrays.copyOf(sj, sj.length);
			}
			else
			{
				numL++;
			}
			if(numL == 500)
				tempValue +=2 ;
			
			numG++;
			
			/* add the new temp value to the yData to plot it in the chart */
//			fitData.add();
//			System.out.println(String.valueOf(fitFunc.apply(getConfiguration(si))));
			bw.write(String.valueOf(fitFunc.apply(getConfiguration(si)))+"\n");
			bw.flush();
		}
		while (numG < maxGlobalIter && numL < maxLokalIter);
		
		bw.close();
		
		System.out.println(String.format("glob : %d, loc : %d", numG,numL));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public int[] findMin() throws SetUpException, InterruptedException
	{
		
		start();
		join();
		return optConf;
	}

	/**
	 * Get a default temperature function
	 */
	public Function<Integer, Double> getDefaultTemperatureFunction()
	{
		return new Function<Integer, Double>()
		{
			@Override
			public Double apply(Integer t)
			{
				return (tempValue *= 0.995);
			}
		};
	}

	/**
	 * Get the found optimized value. This method makes only sense after the
	 * method {@link #findMin()} was run otherwise this method returns null.
	 * 
	 * @return - the optimized value, if the method {@link #findMin()} wasn't
	 *         run null will be returned
	 */
	public R getOptimizedValue()
	{
		if (optConf == null)
			return null;

		return fitFunc.apply(getConfiguration(optConf));
	}

	/**
	 * Get the configuration of the system of the given array of indexes.
	 * 
	 * @param order
	 *            - the array with the indexes
	 * @return - the configuration of the given order of indexes
	 */
	public List<T> getConfiguration(int[] order)
	{
		List<T> currentConf = new ArrayList<>();

		for (int i = 0; i < order.length; i++)
		{
			currentConf.add(objects.get(order[i]));
		}

		return currentConf;
	}

	/**
	 * Generate a random permutation of indexes from 0 to arrayLength-1.
	 * 
	 * @param arrayLength
	 *            - the length of the array (number of elements in the
	 *            permutation)
	 * @return - the random permutation with indexes from 0 to arrayLength and a
	 *         size of arrayLength
	 */
	public int[] generateRandomPermutation(int arrayLength)
	{
		List<Integer> list = new ArrayList<Integer>();
		/*
		 * init the list with natural number of type Integer from 0 to
		 * arrayLength-1
		 */
		for (int i = 0; i < arrayLength; i++)
		{
			list.add(i);
		}

		/* generate random permutation */
		java.util.Collections.shuffle(list);

		/* convert into an array of int */
		int[] permutation = new int[list.size()];
		for (int i = 0; i < permutation.length; i++)
		{
			permutation[i] = list.get(i);
		}

		return permutation;
	}

	/**
	 * Method to generate a new configuration of the system which is "similar"
	 * to the given configuration.
	 * 
	 * @param conf
	 *            - the configuration to change
	 */
	public void changeConfiguration(int[] conf)
	{
		/* init the indexes to change */
		int firstIDX = random(conf.length - 1);
		int secondIDX = random(conf.length - 1);

		/*
		 * generate a new second index until this index is not equal to the
		 * first index
		 */
		while (secondIDX == firstIDX)
		{
			secondIDX = random(conf.length - 1);
		}

		/* switch the positions */
		int temp = conf[firstIDX];
		conf[firstIDX] = conf[secondIDX];
		conf[secondIDX] = temp;

	}

	/**
	 * Generate a natural number between <b>0</b> and <b>max</b>
	 * 
	 * @param max
	 *            - the maximal natural number
	 * @return - an number of type {@code int} between <b>0</b> and <b>max</b>
	 */
	public static int random(int max)
	{
		return (int) (0 + Math.round(Math.random() * max));
	}
}
