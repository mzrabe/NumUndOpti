package org.mzrabe.methaheuristic.firefly;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.None;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;
import org.mzrabe.opti.OptiAlgorithm;
import org.mzrabe.plot.XYContourPlot;

public class FireFlyAlgorithm extends OptiAlgorithm
{
	/**
	 * the list of the flys (agents)
	 */
	private List<Fly> flys = new ArrayList<>();
	
	/**
	 * the function which reduce the intensity over the distance
	 */
	private ReduceIntensity reduceIntensity;
	/**
	 * The n-dimensional range. This ranges will used to generate the random start positions of the Fireflies.
	 * It will also used for plotting.
	 */
	private double[][] range;
	/**
	 * the factor for the random movement
	 */
	private double alpha = 0.05;
	/**
	 * the step length
	 */
	private double beta = 0.6;
	/**
	 * the factor to reduce the influence of the distance
	 */
	private static double gamma = 1e6;
	/**
	 * the number of fireflies (agents or particle)
	 */
	private int numberOfFlies = 150;
	/**
	 * the chart to plot
	 */
	private XYContourPlot chart;
	/**
	 * the directory to plot 
	 */
	private File dir;
	/**
	 * the flag if the steps should plot
	 */
	private boolean shouldPlot = false;
	/**
	 * flag if the tail of the flies should plot
	 */
	private boolean plotTail = true;
	/**
	 * the path where the plot should saved. As default the plots will 
	 * save in the working directory
	 */
	private String pathToSave = "";
	
	/**
	 * last global swarm counter
	 */
	private int lastSwarmCounter = 0;
	/**
	 * formate for the plots
	 */
	private BitmapFormat plotFormat = BitmapFormat.PNG;
	private HashSet<FlySwarm> swarms = new HashSet<>();
	private int numberOfUnchangedIterations = 0;
	private int periode = 15;
	
	/**
	 * @param func
	 *           - the function to optimize (<i>this algorithms search the
	 *            minimum of the given function. If you intent to find a maximum
	 *            you have to multiplied the function by -1</i>)
	 * @param range
	 *            - the range in which the flies will placed randomly
	 * @throws Exception 
	 */
	public FireFlyAlgorithm(Function func, double[][] range) throws Exception
	{
		this(func, getDefaultReduceIntensity(ReduceIntensityFunction.POTENZ), range);
	}
	
	/**
	 * @param func
	 *            - the function to optimize (<i>this algorithms search the
	 *            minimum of the given function. If you intent to find a maximum
	 *            you have to multiplied the function by -1</i>)
	 * @param reduceIntensity
	 *            - the function which describes the reduction of the intensity
	 *            over the distance
	 * @param range
	 *            - the range in which the flies will placed randomly
	 * @throws Exception 
	 */
	public FireFlyAlgorithm(Function func,ReduceIntensity reduceIntensity, double[][] range) throws Exception
	{
		this(func, reduceIntensity, 480, range);
	}
	
	/**
	 * @param func
	 *            - the function to optimize (<i>this algorithms search the
	 *            minimum of the given function. If you intent to find a maximum
	 *            you have to multiplied the function by -1</i>)
	 * @param reduceIntensity
	 *            - the function which describes the reduction of the intensity
	 *            over the distance
	 * @param maxIterations
	 *            - the maximal number of iterations
	 * @param range
	 *            - the range in which the flies will placed randomly
	 * @throws Exception 
	 */
	public FireFlyAlgorithm(Function func,ReduceIntensity reduceIntensity, int maxIterations,double[][] range) throws Exception
	{
		super(func,0.09,maxIterations);
		this.range = range;
		this.reduceIntensity = reduceIntensity;
		
	}
	
	private int getDimension()
	{
		return range.length;
	}

	@Override
	protected void initAlgorithms(double[]... x) throws Exception
	{
		/* distribute the fireflies in the given range */
		for(int i=0;i<numberOfFlies;i++)
		{
			double[] randomPosition = new double[getDimension()];
			
			for(int dim=0;dim<randomPosition.length;dim++)
			{
				double step = Math.abs(range[dim][0]-range[dim][1]);
				randomPosition[dim] = range[dim][0] + step * Math.random();
			}
			
			flys.add(new Fly(randomPosition,i));
			flys.get(i).intensity = -func.getValue(flys.get(i).position);
		}
		
		if (shouldPlot && range.length == 2)
		{
			
			chart = new XYContourPlot("FireFlyAlgorithm - Step 0", "x1", "x2", 800, 600, func, range[0][0], range[0][1], range[1][0], range[1][1], 100,50);
			chart.chart.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);
		    chart.chart.getStyler().setMarkerSize(5);
		    
			dir = new File(pathToSave + "/FireFlyPlots");
			if (!dir.exists())
			{
				dir.mkdirs();
			}

			chart.addDataSet(numberOfFlies + " Flies", getPositions(flys));

			try
			{
				BitmapEncoder.saveBitmap(chart.chart, new File(dir.getAbsolutePath() + String.format("/%1$03d_iteration", numberOfIterations)).getAbsolutePath(), plotFormat);
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
	}

	@Override
	protected void algorithms() throws Exception
	{
		/* make a temporary copy of the current state */
		ArrayList<Fly> copyFly = new ArrayList<>();
//		bestFlys = new ArrayList<>();
		
		for (Fly f : flys)
		{
			copyFly.add(f.copy());
		}
		
		/* unlock the swarms */
		for(FlySwarm swarm : swarms)
		{
			swarm.unlockSwarm();
		}
		
		/** list of the flies which not follows a other fly */
		List<Fly> leaders = new ArrayList<Fly>();

		/* iterate over each fire fly */
		for (int idx_currentFly = 0; idx_currentFly < copyFly.size(); idx_currentFly++)
		{
			/** the current fly to compare with the other */
			Fly currentFly = copyFly.get(idx_currentFly);
			/** the highest intensity, initialize with the intensity with current fly*/
			double bestIntensity = currentFly.intensity;
			/** the fly with the highest intensity */
			Fly bestFly = null;
			/** distance to the other fly which has a higher intensity */
			double distance = 1;
			

			/* find the firefly which have the highest intensity compared with the current fly */
			for (int idx_otherFly = 0; idx_otherFly < copyFly.size(); idx_otherFly++)
			{
				/* don't compare with themselves */
				if(idx_currentFly == idx_otherFly)
					continue;
				
				/** the other fly to compare with the current fly */
				Fly otherFly = flys.get(idx_otherFly);
				/** the geometrical distance between the current fly and the other fly */
				distance = Vector.twoNorm(Vector.minus(otherFly.position, currentFly.position));
				
					/* find the highest intensity and take into account that the intensity decrease over the distance*/
					if ((reduceIntensity.getReducedIntensity(otherFly.intensity, distance) > bestIntensity))
					{
						bestFly = otherFly;
						bestIntensity = reduceIntensity.getReducedIntensity(otherFly.intensity, distance);
					}
			}

			/* generate a random step for the fly */
			double[] random = new double[currentFly.position.length];
			double r = Math.min(Math.abs(range[0][0]-range[0][1]), Math.abs(range[1][0]-range[1][1]));
			for (int idxDim = 0; idxDim < random.length; idxDim++)
			{
				
				random[idxDim] = alpha * (-r/2 + r * Math.random());
			}
			
			/* move the fly */
			if (bestFly == null)
			{
				/* the fly has the highest intensity, so it will only move randomly */
				log.debug(currentFly.toString() + " is the best.");
				currentFly.position = Vector.sum(currentFly.position, random);
				/* save the random step */
				currentFly.randomDirection = random;
				/* add the fly to the leader list */
				leaders.add(currentFly);
			}
			else
			{
				/* the fly sees a other fly with a higher intensity, so it will move a step to the other fly plus a random error */
				log.debug("current fly "+currentFly.ID+" wants to " + bestFly.ID);
				log.debug("current fly "+currentFly.toString());
				currentFly.follows = bestFly.ID;
				
				log.debug("best fly "+bestFly.toString()+"\n");
				
//				/* generate a factor to reduce the random step as closer the fly is on the other fly */
//				double reduceRandom = -Math.pow(1.1,-Math.pow(currentFly.getDistanceToOtherFly(bestFly), 2)) +1.2;
//				random = Vector.multiScalar(random, reduceRandom);
				
				/* generate a factor to simulate the periodic shrinking and expanding of the swarm */
				double period = 15;
				double reduceRandom = Math.sin(numberOfIterations * Math.abs(Math.asin(1.0)/period));
				random = Vector.multiScalar(random, reduceRandom);
				
				
				/* save the step to the other fly */
				currentFly.flowDirection = Vector.multiScalar(Vector.minus(bestFly.position, currentFly.position), beta);
				
				currentFly.position = Vector.sum(currentFly.position, currentFly.flowDirection, random);
				/* save the random step */
				currentFly.randomDirection = random;
				

			
			}
			/* update the intensity of the current fly after taking a step */
			currentFly.intensity = -func.getValue(currentFly.position);
		}
		
		
		HashSet<FlySwarm> currentSwarms = new HashSet<>();
		
		
		/* checks if the leader flies are close to each other */
		for(int i = 0; i < leaders.size();i++)
		{
			/* add the fly to a FlySwarm */
			if(leaders.get(i).follows == null)
			{
				if(leaders.get(i).getFlySwarm() == null || leaders.get(i).getFlySwarm().isLocked())
				{
					/* fly is not in a swarm, so create a new one */
					leaders.get(i).setFlySwarm(new FlySwarm());
				}
	
				
				/* reset the member list of the swarm */
				leaders.get(i).getFlySwarm().members.clear();
				/* update the leaderID of the swarm and a the leader to the member set*/
				leaders.get(i).getFlySwarm().setLeader(leaders.get(i).ID);
				leaders.get(i).getFlySwarm().lockSwarm();
			}
			
			for(int j = i+1; j<leaders.size();j++)
			{
				if(leaders.get(j).follows == null)
				{
					if(leaders.get(i).getDistanceToOtherFly(leaders.get(j)) < 0.5)
					{
						leaders.get(j).follows = leaders.get(i).ID;
					}
				}
			}
		}
		
		/* add all flies to the right swarm */
		for(int i = 0; i<copyFly.size();i++)
		{
			Fly f = copyFly.get(i);
			FlySwarm fsq = getSwarmOfFollowedFly(f, copyFly);
			currentSwarms.add(fsq);
			f.setFlySwarm(fsq);
			

		}
		
		/* plot all swarms */
		int numberFlies = 0;
		for(FlySwarm fsq : currentSwarms)
		{
			log.debug(fsq.toString());
			numberFlies+= fsq.members.size();
			fsq.updatePath(copyFly);
		}
		log.debug("number of flies = " + numberFlies + " in " + currentSwarms.size() + " swarms");
		
		swarms.clear();
		swarms.addAll(currentSwarms);
		
		
		if (shouldPlot)
		{
			/* plot this state as image */
			chart.chart.setTitle("FireFlyAlgorithm - Step " + numberOfIterations);
			chart.updateDataSet(numberOfFlies + " Flies", getPositions(copyFly));
			
			
			
			if (chart.chart.getSeriesMap().get("best flies") == null)
			{
				chart.addDataSet("best flies", getBestPositions(copyFly));
				chart.chart.getSeriesMap().get("best flies").setXYSeriesRenderStyle(XYSeriesRenderStyle.Scatter).setMarkerColor(Color.RED).setLineWidth(4.f).setShowInLegend(false);
			}
			else
			{
				chart.updateDataSet("best flies", getBestPositions(copyFly));
			}
			
			/* plot the swarm path */
			if(true)
			{
				
					for(FlySwarm swarm : currentSwarms)
					{
						List<double[]> path = swarm.path;
						String swarmName = "swarm"+swarm.swarmID;
						if (chart.chart.getSeriesMap().get(swarmName) == null)
						{
						chart.addDataSet(swarmName, path);
						chart.chart.getSeriesMap().get(swarmName).setXYSeriesRenderStyle(XYSeriesRenderStyle.Line).setMarker(new None()).setLineColor(Color.RED).setLineWidth(4.f)
								.setLineStyle(SeriesLines.DASH_DOT).setShowInLegend(false);
						}
						else
						{
							chart.updateDataSet(swarmName, path);
						}
					}
				
			}
			
			if(plotTail)
			{
				/* add the flown direction (tail) */
				for (int number = 0; number < copyFly.size(); number++)
				{
					String flyName_toOther = "fly" + String.format("%1$0" + String.valueOf(copyFly.size()).length() + "d_toOther", number);
					String flyName_random = "fly" + String.format("%1$0" + String.valueOf(copyFly.size()).length() + "d_random", number);
					String flyName_both = "fly" + String.format("%1$0" + String.valueOf(copyFly.size()).length() + "d_both", number);
	
					/** flown path */
					List<double[]> both = new ArrayList<>();
					both.add(flys.get(number).position);
					both.add(copyFly.get(number).position);
	
					/** flown path to other fly */
					List<double[]> toOther = new ArrayList<>();
					toOther.add(copyFly.get(number).getPositionBefor());
					toOther.add(copyFly.get(number).getPositionAfterFlyToOther());
	
					/** random flown path */
					List<double[]> random = new ArrayList<>();
					random.add(copyFly.get(number).getPositionAfterFlyToOther());
					random.add(copyFly.get(number).position);
	
					if (numberOfIterations == 1)
					{
						/* initialize the tails of the fly */
						/* add the flown path to the other fly */
						chart.addDataSet(flyName_both, both);
						/* initialize the style of the path */
						chart.chart.getSeriesMap().get(flyName_both).setXYSeriesRenderStyle(XYSeriesRenderStyle.Line).setMarker(new None()).setLineColor(Color.RED).setLineWidth(1.5f)
								.setLineStyle(SeriesLines.SOLID).setShowInLegend(false);
	
						/* add the random flown path */
						chart.addDataSet(flyName_toOther, toOther);
						/* initialize the style of the path */
						chart.chart.getSeriesMap().get(flyName_toOther).setXYSeriesRenderStyle(XYSeriesRenderStyle.Line).setMarker(new None()).setLineColor(Color.GREEN).setLineWidth(1.0f)
								.setLineStyle(SeriesLines.SOLID).setShowInLegend(false);
	
						/* add the random flown path */
						chart.addDataSet(flyName_random, random);
						/* initialize the style of the path */
						chart.chart.getSeriesMap().get(flyName_random).setXYSeriesRenderStyle(XYSeriesRenderStyle.Line).setMarker(new None()).setLineColor(Color.BLUE).setLineWidth(1.0f)
								.setLineStyle(SeriesLines.SOLID).setShowInLegend(false);
					}
					else
					{
						chart.updateDataSet(flyName_both, both);
						chart.updateDataSet(flyName_toOther, toOther);
						chart.updateDataSet(flyName_random, random);
					}
				}
			}

			try
			{
				BitmapEncoder.saveBitmap(chart.chart, new File(dir.getAbsolutePath() + String.format("/%1$03d_iteration", numberOfIterations)).getAbsolutePath(), plotFormat);
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}

		/* update the current position of every fly */
		for (int j = 0; j < copyFly.size(); j++)
		{
			flys.get(j).intensity = copyFly.get(j).intensity;
			flys.get(j).position = Arrays.copyOf(copyFly.get(j).position, copyFly.get(j).position.length);
			flys.get(j).setFlySwarm(copyFly.get(j).getFlySwarm());
		}
	}
	
	private FlySwarm getSwarmOfFollowedFly(Fly f, List<Fly> flies)
	{
		return f.follows == null ? f.getFlySwarm() : getSwarmOfFollowedFly(flies.get(f.follows), flies);
	}

	@Override
	protected boolean isFinish()
	{
		if(lastSwarmCounter != FlySwarm.nextSwarmID)
		{
			lastSwarmCounter = FlySwarm.nextSwarmID;
//			log.info("nextSwarmID= " + FlySwarm.nextSwarmID);
			numberOfUnchangedIterations = 0;
			return false;
		}
		
		for(FlySwarm swarm : swarms)
		{
			if(swarm.lastSwarmStep() == null || swarm.lastSwarmStep() > tolerance)
			{
//				log.info(numberOfIterations + " - swarm= " + swarm.lastSwarmStep());
				numberOfUnchangedIterations = 0;
				return false;
			}
		}
		
		
		if(numberOfUnchangedIterations++ > 50 && numberOfIterations % (periode*2) == 0)
		{
			return true;
		}
		
		return false;
	}

	@Override
	public String getName()
	{
		return "Fire Fly Algorithm";
	}

	@Override
	public double[] getSolution() throws Exception
	{
		/* find the swarm with the best solution */
		
		double bestValue = Double.MAX_VALUE;
		FlySwarm best = null;
		for(FlySwarm swarm : swarms)
		{
			double value = func.getValue(swarm.getAverage(flys));
			if(value < bestValue)
			{
				bestValue = value;
				best = swarm;
			}
		}
		
		if( best == null)
		{
			/* bilde den schwertpunkt der punkte */
			double[] s = new double[flys.get(0).position.length];
			/* inti with zero */
			Arrays.fill(s, 0.0);
			
			/* calculate the average */
			for (Fly f : flys)
			{
				s = Vector.sum(s, f.position);
			}
			s = Vector.multiScalar(s, 1. / flys.size());
			
			return s;
		}
		else
		{
		return  best.getAverage(flys);
		}
		
	}
	
	public void printAllSolution()
	{
		for(FlySwarm swarm : swarms)
		{
			log.info(swarm.toString());

		}
	}
	
	/**
	 * the interface to define the reduction of the intensity over the distance
	 */
	interface ReduceIntensity
	{
		public double getReducedIntensity(double intensity, double distance);
	}

	/**
	 * @param range the {@link #range} to set
	 */
	public FireFlyAlgorithm setRange(double[][] range)
	{
		int index = 0;
		for(double[] r : range)
		{
			if(r.length % 2 == 0)
			{
				throw new IllegalArgumentException("The range a the index " + index + " is not a valid range. The number values cannot be a odd number.");
			}
			index++;
		}
		
		this.range = range;
		return this;
	}

	/**
	 * The factor the have a influence of the default random step generator.
	 * @param alpha the {@link #alpha} to set
	 */
	public FireFlyAlgorithm setAlpha(double alpha)
	{
		if (alpha > 0)
		{
			this.alpha = alpha;
			return this;
		}
		else
		{
			throw new IllegalArgumentException("The value 'alpha' have to be greater than 0. The value is " + alpha);
		}
	}

	/**
	 * The step length (step in direction of the best firefly). The value has to be in the range ]0,1].
	 * @param beta the {@link #beta} to set
	 */
	public FireFlyAlgorithm setBeta(double beta)
	{
		if (beta > 0 && beta <= 1)
		{
			this.beta = beta;
			return this;
		}
		else
		{
			throw new IllegalArgumentException("The value 'beta' (step length factor) have to be between 0 and 1. The value is " + beta);
		}
	}

	/**
	 * The value must be in the range ]0,1[.
	 * @param gamma the {@link #gamma} to set
	 */
	public FireFlyAlgorithm setGamma(double gamma)
	{
		if (gamma > 0 && gamma < 1)
		{
			FireFlyAlgorithm.gamma = gamma;
			return this;
		}
		else
		{
			throw new IllegalArgumentException("The value 'gamma' have to be between 0 and 1. The value is " + gamma);
		}
	}
	
	/**
	 * Get the positions of the flies as ArrayList
	 */
	public List<double[]> getPositions(List<Fly> flies)
	{
		ArrayList<double[]> list = new ArrayList<>();
		for(Fly f : flies)
		{
			list.add(f.position);
		}
		return list;
	}
	
	/**
	 * Get the positions of the best flies as ArrayList
	 */
	public List<double[]> getBestPositions(List<Fly> flies)
	{
		ArrayList<double[]> list = new ArrayList<>();
		for(Fly f : flies)
		{
			if(f.follows == null)
				list.add(f.position);
		}
		return list;
	}
	
	/**
	 * 
	 * @return - true if the steps will plot
	 */
	public boolean isShouldPlot()
	{
		return shouldPlot;
	}
	/**
	 * Set the flag if the steps should plot or not.
	 * @param shouldPlot - the flag true of false
	 * @return - this for set chains
	 */
	public FireFlyAlgorithm setShouldPlot(boolean shouldPlot)
	{
		this.shouldPlot = shouldPlot;
		return this;
	}
	/**
	 * Get the path where the plots will save
	 * @return - the absolute path to the save directory
	 */
	public String getPathToSave()
	{
		return new File(pathToSave).getAbsolutePath();
	}
	/**
	 * Set the path where the plots will save. In this directory will 
	 * create a directory named 'FireFlyPlots' which contains the images.
	 * This setter checks if the given path is a directory and if it is allowed to write in the directory. 
	 * @param pathToSave - the path (relative or absolute)
	 * @return - this for set chains
	 */
	public FireFlyAlgorithm setPathToSave(String pathToSave)
	{
		File pathToCheck = new File(pathToSave);
		
		if(!pathToCheck.isDirectory())
		{
			throw new IllegalArgumentException("The path "+pathToCheck.getAbsolutePath()+" is not a directory.");
		}
		
		try
		{
			File.createTempFile("test", ".tmp", pathToCheck).deleteOnExit();
			this.pathToSave = pathToSave;
		}
		catch (IOException e)
		{
			throw new IllegalArgumentException("Cannot write at the given path : "+pathToCheck.getAbsolutePath(),e.getCause());
		}
		
		return this;
	}
	
	/**
	 * The Enumerations for the default function approaches which calculate 
	 * the reduction of the intensity depending on the distance.
	 */
	public enum ReduceIntensityFunction
	{
		/**
		 * f(x) = I_0 * exp(-gamma * x^2)
		 */
		EXPONENTIAL,
		/**
		 * f(x) = I_0 / (1 + gamma * x^2)
		 */
		HYPERBEL,
		/**
		 * f(x) = I_0 - x * gamma
		 */
		LINEAR,
		/**
		 * f(x) = I_0
		 */
		NONE, 
		/**
		 * f(x) = I_0 - x^gamma
		 */
		POTENZ
	}
	
	/**
	 * Get the default function to which reduce the intensity over the distance.
	 * @param type - the type of the function approach
	 * @return - the default reduce intensity function
	 */
	private static ReduceIntensity getDefaultReduceIntensity(ReduceIntensityFunction type)
	{
		switch(type)
		{
		case EXPONENTIAL:
			return 	(intensity,distance) -> intensity * Math.exp(-gamma * Math.pow(distance, 1));
		case HYPERBEL:
			return 	(intensity,distance) -> intensity / (1+gamma * Math.pow(distance, 2));
		case LINEAR:
			return 	(intensity,distance) -> -distance*gamma+intensity;
		case POTENZ:
			return  (intensity,distance) -> -Math.pow(distance,gamma)+intensity;
		case NONE:
		default:
			return  (intensity,distance) -> intensity;
		}
	}

	/**
	 * @return the {@link #numberOfFlies}
	 */
	public int getNumberOfFlies()
	{
		return numberOfFlies;
	}

	/**
	 * @param numberOfFlies the {@link #numberOfFlies} to set
	 */
	public void setNumberOfFlies(int numberOfFlies)
	{
		if (numberOfFlies > 1)
		{
			this.numberOfFlies = numberOfFlies;
		}
		else
		{
			throw new IllegalArgumentException("The value 'numberOfFlies' have to be greater than 1. The value is " + numberOfFlies);
		}
	}
}
