
/*	A simulation based on the q-state Potts model
 *	Uses a square lattice with periodic boundary conditions
 *	and the Metropolis algorithm to implement the model
 *  This is a non GUI version for running experiments.
 *  NB temp is in units of Kelvin/Boltzman Constant
 *
 *	Richard Hanes		19/10/04
 */


import java.io.*;

class Lattice
{
    int size, states, equilibriumSweeps, userSweeps, systemEnergy;
    int[] magnetisation, energyArray;
    public int latticeArray[][];
    double temp, finalTemp, tempGradient, heatCapacity, systemMag;
    double magBar, systemMagBar, energyBar, noOfPoints, magSus;
    double[] exponentials, magArray;
    boolean tempChanged = false;
    PrintWriter log;

    public Lattice(int size, int states, double startTemp, double finalTemp, double tempGradient, int sweeps)
    {
	this.size = size;
	this.states = states;
	this.temp = startTemp;
	this.finalTemp = finalTemp;
	this.tempGradient = tempGradient;
	this.userSweeps = sweeps;
	this.equilibriumSweeps = (int)(100 * Math.abs(tempGradient) * ((double)size * (double)size));
	this.latticeArray = new int[size][size];
	this.exponentials = new double[4];
	this.magnetisation = new int[states];
	this.energyArray = new int[sweeps];
	this.magArray = new double[sweeps];
	this.magBar = ((double)(size * size) / (double)states);
	this.noOfPoints = size * size;
	resetSimulation();
	initialiseLog();
    }

    public void resetSimulation()
    {
	resetMagnetArray();
	calculateExponentials();

	if (temp > finalTemp) // If the lattice starts out 'Hot' randomise it
	{                     // otherwise 'Freeze' it by setting it all = 1
	    randomiseLattice();
	} else {
	    orderLattice();
	}
	calculateSystemMagnetisation();
	systemEnergy = calculateTotalSystemEnergy();
    }

    public void resetMagnetArray()
    {
	for (int i = 0; i < states; i++)
	    magnetisation[i] = 0;
    }

    public void calculateExponentials()
    {
	for (int i = 1; i <= 4; i++)
	    exponentials[i-1] = Math.exp(-i/temp);
    }

    public void calculateOverallMagnetisation()
    {
	systemMag = 0;

	for (int i = 0; i < states; i++)
	{
		systemMag += Math.pow((magnetisation[i] - magBar),2);

	   //systemMag += Math.abs(magnetisation[i] - magBar);
	}

	systemMag = Math.sqrt(systemMag);
    }

    public void calculateAverages()
    {
	double eBar = 0; 
	double eSquaredBar = 0;
	double mBar = 0.0;
	double mSquaredBar = 0.0;
	double currentMag = 0.0;
	int currentEnergy = 0;

	for (int i = 0; i < userSweeps; i++){
	    currentEnergy = energyArray[i];
	    eBar += currentEnergy;
	    eSquaredBar += Math.pow(currentEnergy, 2);
	    
	    currentMag = magArray[i];
	    mBar += currentMag;	
	    mSquaredBar += Math.pow(currentMag, 2);		    
	    
	}

	eBar /= userSweeps;
	eSquaredBar /= userSweeps;
	heatCapacity = (eSquaredBar - (eBar * eBar))/ (temp * temp);
	energyBar = eBar;
	
	mBar /= userSweeps;
	mSquaredBar /= userSweeps;
	magSus = (mSquaredBar - (mBar * mBar)) / temp;	
	
	systemMagBar = mBar;

    }


    public void randomiseLattice()
    {
	for (int i = 0; i < size; i++)
	{
	    for (int j = 0; j < size; j++)
	    {
		// Load a random state into each lattice point
		this.latticeArray[i][j] = (int)(Math.random() * states) + 1;
	    }
	}
    }

    public void orderLattice()
    {
	for (int i = 0; i < size; i++)
	{
	    for (int j = 0; j < size; j++)
	    {
		// Load the state 1 into each lattice point
		this.latticeArray[i][j] = 1;
	    }
	}
    }


    public void calculateSystemMagnetisation()
    {
	int currentPointState;

	for (int i = 0; i < size; i++)
	{
	    for (int j = 0; j < size; j++)
	    {
		currentPointState = latticeArray[i][j];
		magnetisation[currentPointState-1]++;
	    }
	}
    }

    // Go through the entire lattice and calculate the energy from each point
    // NB: Do not want to do this each time, waste of processor time
    public int calculateTotalSystemEnergy()
    {
	int energy = 0;
	for (int i = 0; i < size; i++)
	{
	    for (int j = 0; j < size; j++)
	    {
		energy += calculateEnergy(i,j);
	    }
	}
	return energy;
    }

    public void initialiseLog()
    {
	String filename;
	filename = "S";

	filename = filename.concat(String.valueOf(size));
	filename = filename.concat(" Q");
	filename = filename.concat(String.valueOf(states));
	filename = filename.concat(" T");
	filename = filename.concat(String.valueOf(temp));
	filename = filename.concat(" F");
	filename = filename.concat(String.valueOf(finalTemp));
	filename = filename.concat(" G");
	filename = filename.concat(String.valueOf(tempGradient));
	filename = filename.concat(" W");
	filename = filename.concat(String.valueOf(userSweeps));
	filename = filename.concat(".dat");


	try
	{
	    log = new PrintWriter(new FileWriter(filename));

	} catch (Exception e)
	{
	    System.err.println(e.toString());
	    System.exit(1);
	}

	log.println("Temp,Energy,Mag,Cv,MagSus");
    }

    public void outputDataToLog()
    {
	    calculateOverallMagnetisation();
	    calculateAverages();
	
	//log.println(temp + " " + heatCapacity);


	  log.print(temp + "," + energyBar/noOfPoints + "," + systemMagBar/noOfPoints + "," + heatCapacity/noOfPoints + "," + magSus * noOfPoints);
	  for (int i = 0; i < states; i++)
	  {
	  log.print("," + magnetisation[i]/noOfPoints);
	  }
	  log.println("");

    }

    public void outputDataToScreen()
    {
	System.out.println("T: " + temp);// + "\tE: " + energyBar/noOfPoints + "\tCv: " + heatCapacity/noOfPoints + "\tM: " + systemMagBar/noOfPoints);
    }

    // Method to calculate the Energy of a single lattice point
    public int calculateEnergy(int x, int y)
    {
	int[] nearestNeighboursArray = new int[4];
	int targetCellState = latticeArray[x][y];
	int energy = 0;

	//Periodic Boundary Conditions
	if (y == 0)	nearestNeighboursArray[0] = latticeArray[x][size-1];
	else 		nearestNeighboursArray[0] = latticeArray[x][y-1];

	if (x == 0)	nearestNeighboursArray[1] = latticeArray[size-1][y];
	else 		nearestNeighboursArray[1] = latticeArray[x-1][y];

	if (y == (size-1)) 	nearestNeighboursArray[2] = latticeArray[x][0];
	else 				nearestNeighboursArray[2] = latticeArray[x][y+1];

	if (x == (size-1)) 	nearestNeighboursArray[3] = latticeArray[0][y];
	else 				nearestNeighboursArray[3] = latticeArray[x+1][y];

	// This does the work of the Krönecker ð-function
	for (int i = 0; i < 4; i++)
	{
	    if (nearestNeighboursArray[i] == targetCellState)
		energy++;
	}

	return (-1 * energy);
    }

    public void runMetropolis()
    {
	int oldEnergy, newEnergy, energyDifference, systemEnergyDifference;

	// Pick a random lattice point to Flip
	int xToFlip = (int)(Math.random() * size);
	int yToFlip = (int)(Math.random() * size);

	// Remember the old state in case we are gonna flip back
	int oldState = latticeArray[xToFlip][yToFlip];

	// Pick a different random state to flip it to
	int flipToState = (int)(Math.random() * states) + 1;
	while (flipToState == oldState)
	    flipToState = (int)(Math.random() * states) + 1;

	// Calculate old state Energy
	oldEnergy = calculateEnergy(xToFlip,yToFlip);

	// Flip the state
	latticeArray[xToFlip][yToFlip] = flipToState;

	// Calculate the Energy of the new flipped state
	newEnergy = calculateEnergy(xToFlip,yToFlip);
	energyDifference = newEnergy - oldEnergy;

	// The difference in energy of the system is twice the flipped cell difference
	systemEnergyDifference = 2 * energyDifference;

	// If the system is more stable (lower energy) continue
	if (energyDifference < 1)
	{
	    systemEnergy += systemEnergyDifference;
	    magnetisation[oldState-1]--;
	    magnetisation[flipToState-1]++;
	    return;

	} else {
	    // Otherwise flip with a probability based on temperature
	    double flipProb = Math.random();

	    // Only recalculate the exponentials if the temperature has changed
	    if (tempChanged == true)
	    {
		calculateExponentials();
		tempChanged = false;
	    }

	    // Compare the random value flipProb with the Boltzmann Temperature Coefficient
	    if (exponentials[energyDifference-1] > flipProb)
	    {
		systemEnergy += systemEnergyDifference;
		magnetisation[oldState-1]--;
		magnetisation[flipToState-1]++;
		return;

	    } else {
		latticeArray[xToFlip][yToFlip] = oldState;
		return;
	    }
	}
    }
}

public class BlindLatticeSimulator
{
    public static void main (String [] args) throws IOException
    {
	/* Set up Default valuse if none are passed on the Command line */
	int size = 100;
	int states = 5;
	double startTemp = 0.5;
	double finalTemp = 1.5;
	double tempGradient = 0.01;
	int sweeps = 100;
	int i=0;
	int j=0;
	int k=0;
	String argument;

	/* Process Command line Arguments */
	for(int index=0; index < args.length; index++)
	{
	    argument = args[index];
	    if (argument.charAt(0) == '-')
	    {
		char flag = argument.charAt(1);
		argument = argument.substring(2,argument.length());
		switch (flag)
		{
		case 'S':
		    size = Integer.parseInt(argument);
		    break;

		case 'Q':
		    states =  Integer.parseInt(argument);
		    break;

		case 'T':
		    startTemp = Double.parseDouble(argument);
		    break;

		case 'F':
		    finalTemp = Double.parseDouble(argument);
		    break;

		case 'G':
		    tempGradient = Double.parseDouble(argument);
		    break;

		case 'W':
		    sweeps = Integer.parseInt(argument);
		    break;
		}
	    }
	}

	/* Very Basic error checking on initial Conditions */
	if (tempGradient == 0.0)
	{
	    System.out.println("No Temp Gradient, Quitting");
	    System.exit(1);
	}

	if(startTemp < finalTemp)
	{
	    if (tempGradient < 0.0)
	    {
		tempGradient *= -1;
	    }
	}

	if(startTemp > finalTemp)
	{
	    if (tempGradient > 0.0)
	    {
		tempGradient *= -1;
	    }
	}

	/* Initialise the lattice */
	Lattice richie = new Lattice(size, states, startTemp, finalTemp, tempGradient, sweeps);

	/* Print the recieved Initial Conditions to check */
	System.out.println("Input Variables Received:");
	System.out.println("Size:\t" + size);
	System.out.println("States:\t" + states);
	System.out.println("Start Temp:\t" + startTemp);
	System.out.println("Final Temp:\t" + finalTemp);
	System.out.println("Gradient:\t" + tempGradient);
	System.out.println("User Sweeps:\t" + sweeps);
	System.out.println("Equilibrium Sweeps:\t" + richie.equilibriumSweeps);



	if(startTemp < finalTemp)
	{
	    while (richie.temp <= finalTemp)
	    {
		while (i < richie.equilibriumSweeps)
		{
		    richie.runMetropolis();
		    i++;
		}

		if (j >= sweeps)
		{
		    richie.temp += richie.tempGradient;
		    richie.tempChanged = true;
		    richie.outputDataToLog();
		    richie.outputDataToScreen();
		    i = 0;
		    j = 0;
		}

		richie.runMetropolis();
		richie.calculateOverallMagnetisation();
		richie.magArray[j] = richie.systemMag;
		richie.energyArray[j] = richie.systemEnergy;
		j++;
	    }
	    richie.log.close();
	} else {
	    while (richie.temp >= finalTemp)
	    {
	    	while (i < richie.equilibriumSweeps)
		{
		    richie.runMetropolis();
		    i++;
		}

		if (j >= sweeps)
		{
		    richie.temp += richie.tempGradient;
		    richie.tempChanged = true;
		    richie.outputDataToLog();
		    richie.outputDataToScreen();
		    i = 0;
		    j = 0;
		}

		richie.runMetropolis();
		richie.calculateOverallMagnetisation();
		richie.magArray[j] = richie.systemMag;
		richie.energyArray[j] = richie.systemEnergy;
		j++;
	    }
	    richie.log.close();
	}
    }
}

