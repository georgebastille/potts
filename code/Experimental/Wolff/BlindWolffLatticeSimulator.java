
/*
 *  A simulation based on the q-state Potts model
 *  Uses a square lattice with periodic boundary conditions
 *  and the Metropolis algorithm to implement the model
 *  This is a non GUI version for running experiments.
 *  NB temp is in units of Kelvin/Boltzman Constant
 *
 *  Richard Hanes		19/10/04
 */
import java.io.*;

/**
 *  Description of the Class
 *
 * @author     Richie
 * @created    16 November 2004
 */
class Lattice {
	int size, states, equilibriumSweeps, userSweeps, systemEnergy;
	int[] magnetisation, energyArray;
	public int latticeArray[][];
	double temp, finalTemp, tempGradient, addProb, systemMag; 
	double heatCapacity, magSus, magBar, systemMagBar, energyBar, noOfPoints;
	double[] magArray;
	boolean tempChanged = false;
	boolean[][] cluster;
	PrintWriter Log;


	/**
	 *Constructor for the Lattice object
	 *
	 * @param  size          Description of the Parameter
	 * @param  states        Description of the Parameter
	 * @param  startTemp     Description of the Parameter
	 * @param  finalTemp     Description of the Parameter
	 * @param  tempGradient  Description of the Parameter
	 * @param  sweeps        Description of the Parameter
	 */
	public Lattice(int size, int states, double startTemp, double finalTemp, double tempGradient, int sweeps) {
		this.size = size;
		this.states = states;
		this.temp = startTemp;
		this.finalTemp = finalTemp;
		this.tempGradient = tempGradient;
		this.equilibriumSweeps = 10;
		// Maybe calculate this later?
		this.userSweeps = sweeps;
		this.addProb = (1.0 - Math.exp(-1.0 / temp));
		this.latticeArray = new int[size][size];
		this.cluster = new boolean[size][size];
		this.magnetisation = new int[states];
		this.energyArray = new int[sweeps];
		this.magArray = new double[sweeps];
		this.noOfPoints = size * size;
		this.magBar = (noOfPoints / (double) states);

		resetSimulation();
		calculateTotalSystemEnergy();
		calculateSystemMagnetisation();
		initialiseLog();
	}


	/**
	 *  Description of the Method
	 */
	public void resetSimulation() {
		resetMagnetArray();
		resetClusterArray();
		if (temp > finalTemp) {
			// If the lattice starts out 'Hot' randomise it
			// otherwise 'Freeze' it by setting it all = 1
			randomiseLattice();
		} else {
			orderLattice();
		}
	}


	/**
	 *  Description of the Method
	 */
	public void resetMagnetArray() {
		for (int i = 0; i < states; i++) {
			magnetisation[i] = 0;
		}
	}


	/**
	 *  Description of the Method
	 */
	public void calculateOverallMagnetisation() {
		systemMag = 0;

		for (int i = 0; i < states; i++) {

			//systemMag += Math.abs(magnetisation[i] - magBar);
			systemMag += Math.pow((magnetisation[i] - magBar),2);
		}
		systemMag = Math.sqrt(systemMag);

	}


	/**
	 *  Description of the Method
	 */
	public void calculateAverages() {
		double eBar = 0.0;
		double eSquaredBar = 0.0;
		double mBar = 0.0;
		double mSquaredBar = 0.0;
		double currentMag = 0.0;
		int currentEnergy = 0;

		for (int i = 0; i < userSweeps; i++) {
			currentEnergy = energyArray[i];
			eBar += currentEnergy;
			eSquaredBar += Math.pow(currentEnergy, 2);
			
			currentMag = magArray[i];
			mBar += currentMag;
			mSquaredBar += Math.pow(currentMag, 2);
		}

		eBar /= userSweeps;
		eSquaredBar /= userSweeps;
		heatCapacity = (eSquaredBar - (eBar * eBar)) / (temp * temp);
		
		energyBar = eBar;
		
		mBar /= userSweeps;
		mSquaredBar /= userSweeps;
		magSus = (mSquaredBar - (mBar * mBar)) / temp;	
		
		systemMagBar = mBar;
		
	}


	/**
	 *  Description of the Method
	 */
	public void resetClusterArray() {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cluster[i][j] = false;
			}
		}
	}


	/**
	 *  Description of the Method
	 */
	public void randomiseLattice() {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// Load a random state into each lattice point
				this.latticeArray[i][j] = (int) (Math.random() * states) + 1;
			}
		}
	}


	/**
	 *  Description of the Method
	 */
	public void orderLattice() {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// Load the state 1 into each lattice point
				this.latticeArray[i][j] = 1;
			}
		}
	}


	/**
	 *  Description of the Method
	 */
	public void calculateSystemMagnetisation() {
		int currentPointState;

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				currentPointState = latticeArray[i][j];
				magnetisation[currentPointState - 1]++;
			}
		}
	}


	// Go through the entire lattice and calculate the energy from each point
	// NB: Do not want to do this each time, waste of processor time
	/**
	 *  Description of the Method
	 */
	public void calculateTotalSystemEnergy() {
		int energy = 0;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				energy += calculateEnergy(i, j);
			}
		}
		systemEnergy = energy;
	}


	/**
	 *  Description of the Method
	 */
	public void initialiseLog() {
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

		try {
			Log = new PrintWriter(new FileWriter(filename));

		} catch (Exception e) {
			System.err.println(e.toString());
			System.exit(1);
		}

		Log.println("Temp,Energy,Mag,Cv,MagSus");
	}


	/**
	 *  Description of the Method
	 */
	public void calculateObservables() {
		calculateOverallMagnetisation();
		calculateAverages();
	}


	/**
	 *  Description of the Method
	 */
	public void outputDataToLog() {

		Log.print(temp + "," + energyBar / noOfPoints + "," + systemMagBar / noOfPoints + "," + heatCapacity / noOfPoints + "," + magSus * noOfPoints);

		for (int i = 0; i < states; i++) {
			Log.print("," + magnetisation[i] / noOfPoints);
		}
		Log.println("");

	}


	/**
	 *  Description of the Method
	 */
	public void outputDataToScreen() {
		System.out.println("T: " + temp + " E: " + energyBar/ noOfPoints);
	}


	/**
	 *  Description of the Method
	 *
	 * @param  x  Description of the Parameter
	 * @param  y  Description of the Parameter
	 * @return    Description of the Return Value
	 */
	public int calculateEnergy(int x, int y) {
		int[] nearestNeighboursArray = new int[4];
		int targetCellState = latticeArray[x][y];
		int energy = 0;

		//Periodic Boundary Conditions
		if (y == 0) {
			nearestNeighboursArray[0] = latticeArray[x][size - 1];
		} else {
			nearestNeighboursArray[0] = latticeArray[x][y - 1];
		}

		if (x == 0) {
			nearestNeighboursArray[1] = latticeArray[size - 1][y];
		} else {
			nearestNeighboursArray[1] = latticeArray[x - 1][y];
		}

		if (y == (size - 1)) {
			nearestNeighboursArray[2] = latticeArray[x][0];
		} else {
			nearestNeighboursArray[2] = latticeArray[x][y + 1];
		}

		if (x == (size - 1)) {
			nearestNeighboursArray[3] = latticeArray[0][y];
		} else {
			nearestNeighboursArray[3] = latticeArray[x + 1][y];
		}

		// This does the work of the Krönecker ð-function
		for (int i = 0; i < 4; i++) {
			if (nearestNeighboursArray[i] == targetCellState) {
				energy++;
			}
		}

		return (-1 * energy);
	}


	/**
	 *  Description of the Method
	 */
	public void runWolff() {
		if (tempChanged == true) {
			addProb = (1.0 - Math.exp(-1.0 / temp));
			tempChanged = false;
		}

		resetClusterArray();

		int xToFlip = (int) (Math.random() * size);
		int yToFlip = (int) (Math.random() * size);

		int oldState = latticeArray[xToFlip][yToFlip];

		int newState = (int) (Math.random() * states) + 1;
		while (newState == oldState) {
			newState = (int) (Math.random() * states) + 1;
		}

		growCluster(xToFlip, yToFlip, oldState, newState);
	}


	/**
	 *  Description of the Method
	 *
	 * @param  x         Description of the Parameter
	 * @param  y         Description of the Parameter
	 * @param  oldState  Description of the Parameter
	 * @param  newState  Description of the Parameter
	 */
	private void growCluster(int x, int y, int oldState, int newState) {
		cluster[x][y] = true;

		magnetisation[oldState - 1]--;
		magnetisation[newState - 1]++;

		int oldEnergy = calculateEnergy(x, y);

		latticeArray[x][y] = newState;

		int newEnergy = calculateEnergy(x, y);

		int energyDiff = newEnergy - oldEnergy;
		systemEnergy += 2 * energyDiff;

		int xPrev = (x == 0) ? size - 1 : x - 1;
		int xNext = (x == size - 1) ? 0 : x + 1;
		int yPrev = (y == 0) ? size - 1 : y - 1;
		int yNext = (y == size - 1) ? 0 : y + 1;

		if (!cluster[xPrev][y]) {
			tryAdd(xPrev, y, oldState, newState);
		}
		if (!cluster[xNext][y]) {
			tryAdd(xNext, y, oldState, newState);
		}
		if (!cluster[x][yPrev]) {
			tryAdd(x, yPrev, oldState, newState);
		}
		if (!cluster[x][yNext]) {
			tryAdd(x, yNext, oldState, newState);
		}
	}


	/**
	 *  Description of the Method
	 *
	 * @param  x         Description of the Parameter
	 * @param  y         Description of the Parameter
	 * @param  oldState  Description of the Parameter
	 * @param  newState  Description of the Parameter
	 */
	private void tryAdd(int x, int y, int oldState, int newState) {
		if (latticeArray[x][y] == oldState) {
			if (Math.random() < addProb) {
				growCluster(x, y, oldState, newState);
			}
		}
	}
}

/**
 *  Description of the Class
 *
 * @author     Richie
 * @created    16 November 2004
 */
public class BlindWolffLatticeSimulator {


	/**
	 *  The main program for the BlindWolffLatticeSimulator class
	 *
	 * @param  args             The command line arguments
	 * @exception  IOException  Description of the Exception
	 */
	public static void main(String[] args) throws IOException {
		/*
		 *  Set up Default valuse if none are passed on the Command line
		 */
		int size = 32;
		int states = 6;
		double startTemp = 0.5;
		double finalTemp = 5.0;
		double tempGradient = 0.01;
		int sweeps = 1000;
		int i = 0;
		int j = 0;
		int k = 0;
		String argument;

		/*
		 *  Process Command line Arguments
		 */
		for (int index = 0; index < args.length; index++) {
			argument = args[index];
			if (argument.charAt(0) == '-') {
				char flag = argument.charAt(1);
				argument = argument.substring(2, argument.length());
				switch (flag) {
						case 'S':
							size = Integer.parseInt(argument);
							break;
						case 'Q':
							states = Integer.parseInt(argument);
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

		/*
		 *  Initialise the lattice
		 */
		Lattice richie = new Lattice(size, states, startTemp, finalTemp, tempGradient, sweeps);

		/*
		 *  Print the recieved Initial Conditions to check
		 */
		System.out.println("Input Variables Received:");
		System.out.println("Size:\t" + size);
		System.out.println("States:\t" + states);
		System.out.println("Start Temp:\t" + startTemp);
		System.out.println("Final Temp:\t" + finalTemp);
		System.out.println("Gradient:\t" + tempGradient);
		System.out.println("User Sweeps:\t" + sweeps);

		if (startTemp < finalTemp) {
			while (richie.temp <= finalTemp) {

				while (i < richie.equilibriumSweeps) {
					richie.runWolff();
					i++;
				}

				if (j >= sweeps) {
					richie.calculateObservables();
					richie.outputDataToLog();
					richie.outputDataToScreen();
					richie.temp += richie.tempGradient;
					richie.tempChanged = true;
					i = 0;
					j = 0;
				}

				richie.runWolff();
				richie.calculateOverallMagnetisation();
				richie.energyArray[j] = richie.systemEnergy;
				richie.magArray[j] = richie.systemMag;

				j++;
			}
			richie.Log.close();

		} else {
			while (richie.temp >= finalTemp) {
				while (i < richie.equilibriumSweeps) {
					richie.runWolff();
					i++;
				}

				if (j >= sweeps) {
					richie.calculateObservables();
					richie.outputDataToLog();
					richie.outputDataToScreen();
					richie.temp += richie.tempGradient;
					richie.tempChanged = true;
					i = 0;
					j = 0;
				}

				richie.runWolff();
				richie.calculateOverallMagnetisation();
				richie.energyArray[j] = richie.systemEnergy;
				richie.magArray[j] = richie.systemMag;

				j++;
			}
			richie.Log.close();
		}
	}
}

