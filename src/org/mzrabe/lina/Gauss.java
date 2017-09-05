package org.mzrabe.lina;

import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * @author Sebastian Peuser <c-line AT mailbox.tu-berlin DOT de>
 */
public class Gauss {
	
	private static final Logger log = LogManager.getRootLogger();

	/**
	 * Gauss-Jordan-Algorithmus nur fuer eindeutige Gleichungssysteme geeignet
	 * (andernfalls wird NULL zurueckgegeben) matrix[row][column]
	 * @param givenMatrix - the matrix of the linear equation system, coefficient matrix
	 * @param givenVector - the constances of the system
	 * @param printSteps - if the step of the algorithms should print to System.out
	 * @return - the solution vector
	 */
	public static double[] getSolution(double[][] givenMatrix, double[] givenVector, boolean printSteps) {
		
		//make a copy of the matrix
		double[][] matrix = new double[givenMatrix.length][];
		for(int i=0;i<givenMatrix.length;i++)
		{
			matrix[i] = Arrays.copyOf(givenMatrix[i],givenMatrix[i].length);
		}
		//make a copy of the vector
		double[] vector = Arrays.copyOf(givenVector, givenVector.length);
		
		// Das Gleichungssystem hat keine eindeutige Loesung!
		if (matrix.length < matrix[0].length) {
			log.info("Gleichungssystem nicht eindeutig loesbar! The number of rows ("+matrix.length+") are lower than the number of columns ("+matrix[0].length+").");
			return null;
		}

		// Merken der Spalte, welche eine Zahl ungleich null besitzt
		int tmpColumn = -1;

		// Alle Zeilen durchgehen: Ziel der for-Schleife -> Matrix in
		// Zeilenstufenform bringen!
		// -> Alle Zahlen unterhalb der Diagonale sind null
		for (int line = 0; line < matrix.length; line++) {
			tmpColumn = -1;

			// Umformungsschritt 1: Finden einer Spalte mit einem Wert ungleich
			// null
			for (int column = 0; column < matrix[line].length; column++) {
				for (int row = line; row < matrix.length; row++) {
					if (matrix[row][column] != 0) {
						tmpColumn = column;
						break;
					}
				}

				// Abbruch, zahl ungleich null wurde gefunden
				if (tmpColumn != -1) {
					break;
				}
			}

			// NullZeile(n) entdeckt!
			if (tmpColumn == -1) {
				for (int row = line; row < matrix.length; row++) {
					// Gleichungssystem hat keine Loesung!
					if (vector[line] != 0) {
						// Wenn die Zwischenschritte ausgegeben werden sollen
						if (printSteps) {
							printStep(matrix, vector);
						}

						log.info("Gleichungssystem besitzt keine Loesung!");
						return null;
					}
				}
				// Nullzeile(n) vorhanden -> Ist das System noch eindeutig
				// loesbar?
				if (matrix[0].length - 1 >= line) {
					// Wenn die Zwischenschritte ausgegeben werden sollen
					if (printSteps) {
						printStep(matrix, vector);
					}

					// System nicht eindeutig loesbar.
					log.info("Gleichungssystem nicht eindeutig loesbar!");
					return null;
				}
				break;
			}

			// Umformungsschritt 2: Die Zahl matrix[line][tmpColumn] soll
			// UNgleich null sein
			if (matrix[line][tmpColumn] == 0) {
				for (int row = line + 1; row < matrix.length; row++) {
					if (matrix[row][tmpColumn] != 0) {
						// Wenn die Zwischenschritte ausgegeben werden sollen
						if (printSteps) {
							printStep(matrix, vector);
							System.out.println("Zeile " + (line + 1) + " wird mit Zeile " + (row + 1) + " getauscht");
						}

						// Vertauschen von Zeilen -> matrix[line][tmpColumn]
						// wird dann ungleich null
						swapTwoLines(line, row, matrix, vector);
						break;
					}
				}
			}

			// Umformungsschritt 3: matrix[line][tmpColumn] soll gleich 1 sein.
			if (matrix[line][tmpColumn] != 0) {
				// Wenn die Zwischenschritte ausgegeben werden sollen
				if (printSteps) {
					printStep(matrix, vector);
					System.out.println("Zeile " + (line + 1) + " wird durch " + matrix[line][tmpColumn] + " geteilt");
				}

				// Division der Zeile mit matrix[line][tmpColumn]
				divideLine(line, matrix[line][tmpColumn], matrix, vector);
			}

			// Wenn die Zwischenschritte ausgegeben werden sollen
			if (printSteps) {
				printStep(matrix, vector);
			}

			// Umformungsschritt 4: Alle Zahlen unter matrix[line][tmpColumn]
			// sollen null sein.
			for (int row = line + 1; row < matrix.length; row++) {
				// Wenn die Zwischenschritte ausgegeben werden sollen
				if (printSteps) {
					System.out.println("Zu Zeile " + (row + 1) + " wird subtrahiert: " + matrix[row][tmpColumn] + " * Zeile " + (line + 1));
				}

				// Subtraktion damit unter der Zahl im Umformungsschritt 3 nur
				// nullen stehen
				removeRowLeadingNumber(matrix[row][tmpColumn], line, row, matrix, vector);
			}
		}

		// Umformungsschritt 6: Matrix in Normalform bringen (Zahlen oberhalb
		// der Diagonale werden ebenfalls zu null)
		for (int column = matrix[0].length - 1; column > 0; column--) {
			// Wenn die Zwischenschritte ausgegeben werden sollen
			if (printSteps) {
				printStep(matrix, vector);
			}

			// Alle Werte oberhalb von "column" werden zu null
			for (int row = column; row > 0; row--) {
				// Wenn die Zwischenschritte ausgegeben werden sollen
				if (printSteps) {
					System.out.println("Zu Zeile " + (row) + " wird subtrahiert: " + matrix[row - 1][column] + " * Zeile " + (column + 1));
				}

				// Dazu wird Subtraktion angewandt
				removeRowLeadingNumber(matrix[row - 1][column], column, row - 1, matrix, vector);
			}
		}

		// Wenn die Zwischenschritte ausgegeben werden sollen
		if (printSteps) {
			printStep(matrix, vector);
		}

		// Unser ehemaliger Loesungsvektor ist jetzt zu unserem Zielvektor
		// geworden :)
		return vector;
	}

	/**
	 * Hier werden einfach zwei Zeilen vertrauscht
	 */
	private static void swapTwoLines(int rowOne, int rowTwo, double[][] matrix, double[] vector) {
		double[] tmpLine;
		double tmpVar;

		tmpLine = matrix[rowOne];
		tmpVar = vector[rowOne];

		matrix[rowOne] = matrix[rowTwo];
		vector[rowOne] = vector[rowTwo];

		matrix[rowTwo] = tmpLine;
		vector[rowTwo] = tmpVar;
	}

	/**
	 * eine Zeile wird durch "div" geteilt. "div" darf nicht null sein
	 */
	private static void divideLine(int row, double div, double[][] matrix, double[] vector) {
		for (int column = 0; column < matrix[row].length; column++) {
			matrix[row][column] = matrix[row][column] / div;
		}
		vector[row] = vector[row] / div;
	}

	/**
	 * Eine Zeile (row) wird mit einem entsprechendem vielfachen (factor) von
	 * einer anderen Zeile (rowRoot) subtrahiert.
	 */
	private static void removeRowLeadingNumber(double factor, int rowRoot, int row, double[][] matrix, double[] vector) {
		for (int column = 0; column < matrix[row].length; column++) {
			matrix[row][column] = matrix[row][column] - factor * matrix[rowRoot][column];
		}
		vector[row] = vector[row] - factor * vector[rowRoot];
	}

	/**
	 * Ein Vector wird auf der Konsole ausgegeben (transponiert)
	 * @param vector - the vector to print
	 */
	public static void printVector(double[] vector) {
		
		if (vector == null) {
			return;
		}
		System.out.print("Loesungsvektor ist: (");
		for (int i = 0; i < vector.length; i++) {
			if (i != 0) {
				System.out.print(",");
			}
			System.out.print(vector[i]);
		}
		System.out.println(")^T");
	}

	/**
	 * Eine Matrix wird auf der Konsole ausgegeben matrix[row][column]
	 * @param matrix - the matrix to print
	 */
	public void printMatrix(double[][] matrix) {
		if (matrix == null) {
			return;
		}
		for (int row = 0; row < matrix.length; row++) {
			System.out.print("(");
			for (int column = 0; column < matrix[row].length; column++) {
				if (column != 0) {
					System.out.print(",");
				}
				System.out.print(matrix[row][column]);
			}
			System.out.println(")");
		}
	}

	/**
	 * Diese Methode zeigt die Zwischenschritte der Berechnung auf der Konsole
	 * an. Fuer die Aufgabe nicht weiter relevant (unbekannte Konzepte werden
	 * verwendet!)
	 */
	private static void printStep(double[][] matrix, double[] vector) {
		System.out.println();

		// Werte werden fuer die Ausgabe auf ein bestimmtes Format gebracht
		// -> Damit die Ausgabe auch immer schick aussieht
		java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");
		for (int row = 0; row < matrix.length; row++) {
			for (int column = 0; column < matrix[row].length; column++) {
				if (matrix[row][column] >= 0) {
					System.out.print("+");
				}
				System.out.print(df.format(matrix[row][column]) + "    ");
			}
			System.out.print("|    ");
			if (vector[row] >= 0) {
				System.out.print("+");
			}
			System.out.println(df.format(vector[row]));
		}
	}
	
}
