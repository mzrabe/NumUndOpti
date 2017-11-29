package org.mzrabe.opti;

import static org.mzrabe.lina.Vector.multiScalar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.MathOperation;
import org.mzrabe.lina.Vector;

/**
 * @author Moritz Zahn, email: mzrabe@gmail.com <br><br>
 *
 */
public class PenaltyMethod extends OptiAlgorithm
{
	
	public interface IncreaseDefinition
	{
		/**
		 * The mathematical operation to increase the penalty coefficient after violent the restriction.
		 * @param penaltyCoefficient - the penaltyCoefficient to increase
		 */
		public double increase(double penaltyCoefficient);
	}
	
	/**
	 * The constraint for a optimization problem. <br>
	 * Inequalities must have the form g(x) <= 0 and equation 
	 * must have the form h(x) = 0 where x is an n dimensional vector.
	 */
	public class Restriction
	{
		/**
		 * the type of the restriction
		 */
		public RestrictionsType type;
		/**
		 * The definition of the restriction. Inequalities must have the form g(x) <= 0 and equation 
		 * must have the form h(x) = 0.
		 */
		public MathOperation restrictionDefinition;
		/**
		 * The mathematical operation to increase the penalty coefficient after violent the restriction.
		 */
		public IncreaseDefinition increaseDefinition;
		/**
		 * the penalty coefficient
		 */
		public double r;
		/**
		 * is true if the restriction was violent
		 */
		public boolean violent;
		
		/**
		 * Constructor for a restriction for the penalty optimization method.
		 * 
		 * @param type - the type of the restriction
		 * @param r - the start penalty coefficient
		 */
		public Restriction(RestrictionsType type, double r)
		{
			this.type = type;
			this.r = r;
			this.violent = false;
		}
		
		/**
		 * @param x - the function arguments of the defined restriction function
		 * @param c - a optional number of changeable coefficients of the defined restriction function
		 * @return the function value of the restriction
		 */
		public Double getPenaltyValue(double[] x, double ... c )
		{
			switch(this.type)
			{
				case EQUATION:
					violent = true;
					return Math.pow(restrictionDefinition.getValue(x, c), 2) * r;
				case INEQUALITY:
					double value = restrictionDefinition.getValue(x, c);
//					log.trace("x = " + Arrays.toString(x) + " value = " + value);
					if(value > 0)
					{
						violent = true;
						return Math.pow(restrictionDefinition.getValue(x, c), 2) * r;
					}
					else
					{
						return 0.0;
					}
				default:
					System.out.println("The type of the restriction was not set. Null is returned.");
					return null;
			}
		}
		
		/**
		 * Increase the penalty coefficient under respect of the 
		 * {@link RestrictionDefinition#increaseDefinition(double)} it the restriction was violent.
		 */
		public void increasePenaltyCoefficient()
		{
			if(violent)
			{
				r = increaseDefinition.increase(r);
				violent = false;
			}
		}
		
		/**
		 * set the definition of the restriction
		 */
		public Restriction setRestrictionDefinition(MathOperation restrictionDefinition)
		{
			this.restrictionDefinition = restrictionDefinition;
			return this;
		}
		/**
		 * set the definition of the mathematical operation to increase the penalty coefficient
		 */
		public Restriction setIncreaseDefinition(IncreaseDefinition increaseDefinition)
		{
			this.increaseDefinition = increaseDefinition;
			return this;
		}
	}
	
	/**
	 * the types of the restrictions
	 *
	 */
	public enum RestrictionsType
	{
		INEQUALITY,EQUATION
	}

	/**
	 * The list of the restriction like equation or inequalities.
	 */
	private List<Restriction> restrictions = new ArrayList<>();
	/**
	 * The original function without the restrictions (constraints)
	 */
	@SuppressWarnings("unused")
	private Function originFunction;
//	/**
//	 * The penalty function which include the constraints.
//	 */
//	public Function penaltyFunc;
	/**
	 * the sub optimization algorithms
	 */
	private OptiAlgorithm optiAlgo;
	/**
	 * the last minimum
	 */
	private double[] lastMin;
	/**
	 * the start points of the simplex
	 */
	private double[][] startSimplex;

	/**
	 * Constructor to init the penalty method
	 * 
	 * @param restrictions - a list of the restriction for the optimization problem
	 * @param originFunction - the function which describes the optimization problem
	 * @param tolerance - the tolerance for the algorithms
	 * @param maxIterations - the maximal iterations for the algorithms
	 * @param optiAlgo - the optimization algorithms to find the minimum of the penalty function inside the iterations
	 */
	public PenaltyMethod(List<Restriction> restrictions, Function originFunction, double tolerance, int maxIterations, OptiAlgorithm optiAlgo)
	{
		super(tolerance, maxIterations);
		this.restrictions = restrictions;
		this.optiAlgo = optiAlgo;
		initPenaltyFunc(originFunction);
		
	}
	
	/**
	 * Constructor to init the penalty method. If this constructor is used the restrictions has to be set
	 * 
	 * @param originFunction - the function which describes the optimization problem
	 * @param tolerance - the tolerance for the algorithms
	 * @param maxIterations - the maximal iterations for the algorithms
	 * @param optiAlgo - the optimization algorithms to find the minimum of the penalty function inside the iterations
	 */
	public PenaltyMethod(Function originFunction, double tolerance, int maxIterations, OptiAlgorithm optiAlgo)
	{
		super(tolerance, maxIterations);
		this.optiAlgo = optiAlgo;
		initPenaltyFunc(originFunction);
		
	}
	
	
	/**
	 * Initialized the penalty function.
	 * @param originFunction - the original function to optimize
	 * @return this instance
	 */
	public PenaltyMethod initPenaltyFunc(Function originFunction)
	{
		
		this.originFunction = originFunction;
		this.func = new Function()
		{

			@Override
			public double getValue(double[] x, double... c)
			{

				double value = originFunction.getValue(x, c);
				if(restrictions.size() == 0)
				{
					log.warn("No restrictions are set.");
				}
				for(Restriction res : restrictions)
				{
//					log.trace("res : " +i);
					value += res.getPenaltyValue(x, c);
				}

				return value;
			}
		};
		
		this.optiAlgo.setFunc(func);
		return this;
	}

	@Override
	protected void initAlgorithms(double[]... x)
	{
//		if(optiAlgo instanceof DownhillSimplex)
//		{
//			lastMin = ((DownhillSimplex) optiAlgo).getPsi_();
//		}
//		else
//		{
			lastMin = Arrays.copyOf(x[0], x[0].length);
//	}
		startSimplex = x;
		addPointsToHistory(lastMin);
		
		log.debug("run with penalty coefficient(s) of : " + getPCO());
		
		try
		{
			solution = optiAlgo.findMin(startSimplex);
		}
		catch (InterruptedException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		increasePenaltyCoefficients();
	}

	@Override
	protected void algorithms()
	{
		lastMin = Arrays.copyOf(solution, solution.length);
		
		log.debug("run with penalty coefficient(s) of : " + getPCO());
		if(optiAlgo instanceof DownhillSimplex)
		{
			try
			{
				solution = optiAlgo.findMin(startSimplex);
			}
			catch (InterruptedException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else
		{
			try
			{
				solution = optiAlgo.findMin(lastMin);
			}
			catch (InterruptedException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
		}
		
		increasePenaltyCoefficients();
		
	}

	@Override
	protected boolean isFinish()
	{
//		log.debug("solution: " + Arrays.toString(solution) + "lastMin: " + Arrays.toString(lastMin) + " tol: " +Vector.oneNorm(Vector.minus(solution, lastMin)));
		
		return Vector.oneNorm(Vector.minus(solution, lastMin)) < tolerance ? true : false;
	}

	@Override
	public String getName()
	{
		return "Penalty Method (with " + (restrictions != null ? restrictions.size() : 0) + " restrictions)";
	}
	
	/**
	 * Increase the penalty coefficients. 
	 * This method has only a effect on the penalty coefficient if the restriction was violent. 
	 */
	protected void increasePenaltyCoefficients()
	{
		for(Restriction res : restrictions)
		{
			res.increasePenaltyCoefficient();
		}
	}
	
	/**
	 * Add a restriction to the penalty optimization method.
	 * @param restriction - a restriction like a equation or inequality of type {@code Restriction}
	 */
	public PenaltyMethod addRestriction(Restriction restriction)
	{
		restrictions.add(restriction);
		return this;
	}
	
	/**
	 * 
	 * @param type - the type of the restriction
	 * @param penaltyCoefficient - the start value of the penalty coefficient
	 * @param restrictionDefinition - the definition of the restriction
	 * @param increaseDefinition - the definition to increase the penalty coefficient
	 */
	public PenaltyMethod addRestriction(RestrictionsType type, double penaltyCoefficient, MathOperation restrictionDefinition, IncreaseDefinition increaseDefinition)
	{
		restrictions.add(new Restriction(type, penaltyCoefficient).setRestrictionDefinition(restrictionDefinition).setIncreaseDefinition(increaseDefinition));
		return this;
	}
	
	private String getPCO()
	{
		String str = "[ ";
		for(Restriction r : restrictions)
		{
			str+=r.r + " ";
		}
		
		str+= "]";
		return str;
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}
	
	
	
	

}
