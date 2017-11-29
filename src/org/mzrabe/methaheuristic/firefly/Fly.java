package org.mzrabe.methaheuristic.firefly;

import java.util.Arrays;
import java.util.HashSet;

import org.mzrabe.lina.Vector;

/**
 *
 */
public class Fly
{
	/**
	 * the id of the fly
	 */
	public int ID;
	/**
	 * the intensity of the fly
	 */
	public double intensity;
	/**
	 * the position of the fly in the n-dimensional space
	 */
	public double[] position;
	/**
	 * the flow direction (vector) to the other (better) fly
	 */
	public double[] flowDirection = null;
	/**
	 * the random direction (vector) of the fly
	 */
	public double[] randomDirection = null;
	/**
	 * the swarm in which the fly is a member (the group)
	 */
	private FlySwarm flySwarm;
//	/**
//	 * the flies which follows this fly, means that this fly is the most attractive for the other flies in this set
//	 */
//	public HashSet<Fly> followers = new HashSet<>();
	public Integer follows = null;
	
	/**
	 * Constructor to get a new instance of a fly.
	 * @param startPosition - the start position of the fly
	 * @param id - the id of the fly
	 */
	public Fly(double[] startPosition,int id)
	{
		position = startPosition;
		this.ID = id;
	}
	
	/**
	 * Empyt Constructor
	 */
	public Fly(){}
	
	/**
	 * Get a deep copy of the fly.
	 * @return - a deep copy of this fly
	 */
	public Fly copy()
	{
		Fly copy = new Fly();
		copy.ID = this.ID;
		copy.intensity = this.intensity;
		copy.position = Arrays.copyOf(this.position, this.position.length);
		copy.flySwarm = this.flySwarm;
		if(flowDirection != null)
			copy.flowDirection = Arrays.copyOf(this.flowDirection, this.flowDirection.length);
		if(randomDirection != null)
			copy.randomDirection = Arrays.copyOf(this.randomDirection, this.randomDirection.length);
		return copy;
	}
	
	/**
	 * Get the last position of the fly.
	 * @return - the position before the current position
	 */
	public double[] getPositionBefor()
	{
		double[] pos = Arrays.copyOf(position, position.length);
		
		if(flowDirection != null)
			pos = Vector.minus(pos, flowDirection);
		if(randomDirection != null)
			pos = Vector.minus(pos, randomDirection);
		
		return pos;
	}
	/**
	 * Get the position before the fly made a random step. Which mean the current position minus the random step.
	 * @return - the position of the fly before random step
	 */
	public double[] getPositionAfterRandom()
	{
		double[] pos = Arrays.copyOf(getPositionBefor(), position.length);
		
		if(randomDirection != null)
			pos = Vector.sum(pos, randomDirection);
		
		return pos;
	}
	/**
	 * Get the position before the fly made a step to the better (more attractive) fly. The current position minus the flown step to the other fly.
	 * @return - the position before made a step the other fly
	 */
	public double[] getPositionAfterFlyToOther()
	{
		double[] pos = Arrays.copyOf(getPositionBefor(), position.length);
		
		if(flowDirection != null)
			pos = Vector.sum(pos, flowDirection);
		
		return pos;
	}
	/**
	 * @return - true if the fly made a step to a other (better/ more attractive) fly. Otherwise false;
	 */
	public boolean flowToOther()
	{
		return flowDirection == null ? false : true;
	}
	/**
	 * @return -true if the fly made a random step
	 */
	public boolean flowRandom()
	{
		return randomDirection == null ? false : true;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
//		Integer[] fIds = new Integer[followers.size()];
//		int i = 0;
//		for(Fly f : followers)
//		{
//			fIds[i] = f.ID;
//			i++;
//		}
		
				
		return String.format("id= %d, intensity= %f, position= %s, flowDirection= %s, randomDirection= %s, follows= %s", ID, intensity,Arrays.toString(position),Arrays.toString(flowDirection),Arrays.toString(flowDirection), follows == null ? "NOBODY" : follows.toString());
	}
	
	/**
	 * Get the distance (Euclidean norm) from this fly to the other (given) fly.
	 * @param otherFly - the other fly 
	 * @return - the distance between this and the other fly
	 */
	public double getDistanceToOtherFly(Fly otherFly)
	{
		return Vector.twoNorm(Vector.minus(position, otherFly.position));
	}

	
	/**
	 * Get the swarm in which the fly is a member.
	 * @return - the swarm in which the fly is a member
	 */
	public FlySwarm getFlySwarm()
	{
		return this.flySwarm;
	}
	
	/**
	 * Set the swarm in which the fly is a member.
	 * @param - the swarm to set
	 */
	public void setFlySwarm(FlySwarm flySwarm)
	{
		this.flySwarm = flySwarm;
		this.flySwarm.members.add(this.ID);
	}
	
	
}
