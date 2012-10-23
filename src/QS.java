import java.math.BigInteger;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

public class QS {

	// Constants
	private final BigInteger ZERO = BigInteger.ZERO;
	private final BigInteger ONE = BigInteger.ONE;
	private final BigInteger TWO = ONE.add(ONE);
	private int iterationLimit;
	private int flattenFactor;
	private int PRIME_BASE;

	// Fields
	private BigInteger[] bigPrimes;
	private int[] primes;
	private double[] primeLogs;
	private BigInteger N;
	BigInteger[] s;
	BigInteger[] fs;
	double[] logs;
	Vector decomposedNumbers;
	int decomposedCount;

	public static void main(String[] args) {
		for (String s : args) {
			List<BigInteger> list = getFactors(new BigInteger(s));
			for (BigInteger bi : list) {
				System.out.println(bi);
			}
		}
	}

	public static List<BigInteger> getFactors(BigInteger N) {
		List<BigInteger> factors = new LinkedList<BigInteger>();
		while (!N.isProbablePrime(50)) {
			QS qs = new QS(N);
			BigInteger factor = qs.quadraticSieve(0);
			if (factor.isProbablePrime(50)) {
				factors.add(factor);
				N = N.divide(factor);
			} else {
				List<BigInteger> facts = getFactors(factor);
				factors.addAll(getFactors(factor));
				for (BigInteger bi : facts) {
					N = N.divide(bi);
				}
			}
		}
		factors.add(N);
		return factors;
	}

	private static void println(String s) {
		System.out.println(s);
	}

	public QS(BigInteger toFactor) {
		this.N = toFactor;

		// flattenFactor = 20;

		// double logn = N.bitLength() * Math.log(2);
		// PRIME_BASE = (int) (Math.pow(Math.E, 0.5 * Math.sqrt((logn *
		// Math.log(logn)))) / 2);
		//
		// if (PRIME_BASE > 1500)
		// PRIME_BASE = 1500;
		//
		// while (PRIME_BASE % flattenFactor != 0) {
		// PRIME_BASE++;
		// }
		// iterationLimit = 20000;
		int bitCount = N.bitLength();
		if (bitCount < 85) {
			PRIME_BASE = 500;
			flattenFactor = 10;
			iterationLimit = 15000;
		} else {
			PRIME_BASE = 1000;
			flattenFactor = 20;
			iterationLimit = 20000;
		}

		// Initializes an array of the first <PRIME_BASE> prime numbers, in
		// BigInteger form.
		bigPrimes = new BigInteger[PRIME_BASE];
		primeLogs = new double[PRIME_BASE];
		primes = new int[PRIME_BASE];
		primes[0] = -1;
		bigPrimes[0] = BigInteger.valueOf(-1);
		primeLogs[0] = 0;
		int i = 1;
		int j = 1;
		while (j < PRIME_BASE) {
			// Analyzes the next prime from the raw array.
			int testPrime = Primes_pollard.primes[i];
			BigInteger bigPrime = BigInteger.valueOf(testPrime);

			// Calculates N % p, where p is our prime for testing.
			int nModP = N.mod(bigPrime).intValue();

			// Checks if the prime is a Quadratic residue of N.
			if (legendreSymbol(nModP, testPrime) == 1) {
				primes[j] = testPrime;
				bigPrimes[j] = bigPrime;
				primeLogs[j] = Math.log(testPrime);
				j++;
			}
			i++;
		}
	}

	public BigInteger quadraticSieve(long looplimit) {
		// The size of various computational vectors.

		// Smooth number pairs
		int smoothNumbers[][] = new int[PRIME_BASE][2];

		// Vector of decomposed numbers.
		decomposedNumbers = new Vector();
		decomposedCount = 0;

		// Congruences
		for (int i = 1; i < PRIME_BASE; i++) {
			// Finds smooth numbers over the factor base.
			long[] smooth = findQxPolynomials(primes[i], N.mod(bigPrimes[i]).intValue());
			smoothNumbers[i][0] = (int) smooth[0];
			smoothNumbers[i][1] = (int) smooth[1];
		}

		// m = The root of N, rounded up.
		BigInteger m = sqrt(N);

		// This will be a large vectors of "random" numbers
		s = new BigInteger[iterationLimit + 2];

		// This will be the "random" numbers squared modulo
		fs = new BigInteger[iterationLimit + 2];

		logs = new double[iterationLimit + 2];

		int direction = 0;
		int offset = 0;
		do {
			// Decides parameters for finding the next s/fs-values.
			switch (direction) {
			case 0:
				direction = 1;
				offset = 0;
				break;
			case 1:
				direction = -1;
				offset = -offset + iterationLimit;
				break;
			case -1:
				direction = 1;
				offset = -offset;
				break;
			}

			// Cleans the s/vs/logs-arrays.
			for (int i = 0; i < iterationLimit; i++) {
				s[i] = null;
				fs[i] = null;
				logs[i] = 0;
			}

			// Recalculates the s, sf and log vectors.
			for (int i = 1; i < PRIME_BASE; i++) {

				// The prime to look at.
				int p = primes[i];

				// The logarithm of this prime.
				double logp = primeLogs[i];

				// The M-value % p + offset.
				int mInt = m.mod(bigPrimes[i]).intValue() + offset;

				// Checks if the congruence is solved.
				if (smoothNumbers[i][0] >= 0) {
					int loop = (int) ((smoothNumbers[i][0] - mInt) % p);
					if (loop < 0)
						loop += p;
					// Sieves in the logs-array.
					while (loop < iterationLimit) {
						logs[loop] += logp;
						loop += p;
					}
				}
				// Checks if the congruence is solved.
				if (smoothNumbers[i][1] >= 0) {
					int loop = (int) ((smoothNumbers[i][1] - mInt) % p);
					if (loop < 0)
						loop += p;
					// Sieves in the logs-array.
					while (loop < iterationLimit) {
						logs[loop] += logp;
						loop += p;
					}
				}
			}

			// The value that a sieved element must have to have been hit enough
			// times.
			double TARGET = (double) (Math.log(m.doubleValue()) + Math.log(iterationLimit) - primeLogs[PRIME_BASE - 1]);

			for (int i = 0; i < iterationLimit; i++) {
				// Checks the sieve for target value (if it has been hit by the
				// sieve enough times).
				if (logs[i] > TARGET) {
					// Calculates s[i] and fs[i] if the sieve hit this number.
					s[i] = BigInteger.valueOf(i + offset).add(m);
					fs[i] = s[i].multiply(s[i]).subtract(N);
				}
			}

			for (int i = 1; i < PRIME_BASE; i++) {

				int prime = primes[i];
				int mInt = m.mod(bigPrimes[i]).intValue() + offset;

				// Checks if the congruence is solved.
				if (smoothNumbers[i][0] >= 0) {
					int loop = (int) ((smoothNumbers[i][0] - mInt) % prime);
					if (loop < 0)
						loop += prime;
					while (loop < iterationLimit) {
						if (logs[loop] > TARGET) {
							divideByPower(loop, bigPrimes[i]);
							if (decomposedCount >= PRIME_BASE + flattenFactor)
								break;
						}
						loop += prime;
					}
				}
				if (smoothNumbers[i][1] >= 0) {
					int loop = (int) ((smoothNumbers[i][1] - mInt) % prime);
					if (loop < 0)
						loop += prime;
					// Going through the sieve to decompose numbers.
					while (loop < iterationLimit) {
						if (logs[loop] > TARGET) {
							// Divides f[loop] with the current prime as long as
							// it's possible. If it ends up with 1, it is fully
							// decomposed and added to the vector.
							divideByPower(loop, bigPrimes[i]);
							// Enough numbers are decomposed, ends the looping.
							if (decomposedCount >= PRIME_BASE + flattenFactor)
								break;
						}
						loop += prime;
					}
				}
				if (decomposedCount >= PRIME_BASE + flattenFactor)
					break;
			}

			// Keeps looping until enough congruences are checked and
			// decomposed.
		} while (decomposedCount < PRIME_BASE + flattenFactor);

		// If too many decomposed numbers are found, the topmost ones will be
		// ignored.
		if (decomposedCount > PRIME_BASE + flattenFactor)
			decomposedCount = PRIME_BASE + flattenFactor;

		// Recalculates s and f(s) according to our decomposed numbers.
		s = new BigInteger[decomposedCount + 1];
		fs = new BigInteger[decomposedCount + 1];
		for (int i = 0; i < decomposedCount; i++) {
			s[i] = ((BigInteger) decomposedNumbers.elementAt(i));
			fs[i] = s[i].multiply(s[i]).subtract(N);
		}
		decomposedNumbers = null;

		// Creates the identity matrix (flattened)
		int id[][] = buildIdentityMatrix(decomposedCount);
		// Creates the coef matrix
		byte factors[][] = buildMatrix(fs, decomposedCount);
		// Flattens the coef matrix to match the identity.
		int matrix[][] = flattenMatrix(factors);

		// Solving the linear equation by gauss-elimination.
		solveMatrix(matrix, id);
		BigInteger test = ONE;

		// Loops through the decomposed vector from the top down, to find rows
		// which correspond to quadratic congruences.
		for (int loop = decomposedCount - 1; loop > PRIME_BASE; loop--) {
			int primeFactors[] = new int[PRIME_BASE];
			byte factorLine[] = matrixLine(id, loop);

			test = ONE;
			// Initializes the prime factors to 0.
			for (int i = 0; i < PRIME_BASE; i++)
				primeFactors[i] = 0;
			for (int i = 0; i < decomposedCount; i++) {
				if (factorLine[i] == 1) {
					for (int j = 0; j < PRIME_BASE; j++)
						primeFactors[j] += (int) factors[i][j];
					test = test.multiply(s[i]).mod(N);
				}
			}

			BigInteger prim = ONE;
			for (int i = 0; i < PRIME_BASE; i++) {
				BigInteger y = BigInteger.valueOf(primes[i]).modPow(BigInteger.valueOf(primeFactors[i] / 2), N);
				prim = prim.multiply(y).mod(N);
			}

			test = test.mod(N);
			prim = prim.mod(N);

			BigInteger x = test.add(prim);
			BigInteger y = test.subtract(prim);

			test = N.gcd(x);
			if ((test.compareTo(ONE) != 0) && (test.compareTo(N) != 0))
				break;
		}

		return test;
	}

	private byte[] matrixLine(int matrix[][], int index) {
		int j, k;
		int line[] = matrix[index];
		byte result[] = new byte[PRIME_BASE + flattenFactor];
		int comparation = 1;
		for (j = 0; j < PRIME_BASE + flattenFactor; j++) {
			if ((line[j / flattenFactor] & comparation) > 0)
				result[j] = 1;
			else
				result[j] = 0;

			if (j % flattenFactor == (flattenFactor - 1))
				comparation = 1;
			else
				comparation *= 2;
		}
		return result;
	}

	private void gaussElim(int matrix[][], int right[][], int j, int k) {
		int c1, c2;
		int temp;
		int comparation = 1;
		for (c1 = 1; c1 <= (j % flattenFactor); c1++)
			comparation *= 2;
		if ((matrix[j][j / flattenFactor] & comparation) == 0) {
			for (c1 = j + 1; c1 < k; c1++) {
				if ((matrix[c1][j / flattenFactor] & comparation) > 0) {
					for (c2 = j / flattenFactor; c2 < k / flattenFactor; c2++) {
						temp = matrix[j][c2];
						matrix[j][c2] = matrix[c1][c2];
						matrix[c1][c2] = temp;
					}
					for (c2 = 0; c2 < k / flattenFactor; c2++) {
						temp = right[j][c2];
						right[j][c2] = right[c1][c2];
						right[c1][c2] = temp;
					}
					break;
				}
			}
		}

		if ((matrix[j][j / flattenFactor] & comparation) > 0) {
			for (c1 = j + 1; c1 < k; c1++) {
				if ((matrix[c1][j / flattenFactor] & comparation) > 0) {
					for (c2 = j / flattenFactor; c2 < k / flattenFactor; c2++)
						matrix[c1][c2] = (int) (matrix[c1][c2] ^ matrix[j][c2]);
					for (c2 = 0; c2 < k / flattenFactor; c2++)
						right[c1][c2] = (int) (right[c1][c2] ^ right[j][c2]);
				}
			}
		}
	}

	private void solveMatrix(int matrix[][], int right[][]) {
		int j, k;
		k = matrix.length;
		for (j = 0; j < k - 1; j++)
			gaussElim(matrix, right, j, k);
	}

	private byte[][] buildMatrix(BigInteger numbers[], int size) {
		byte matrix[][] = new byte[size][size];
		BigInteger temp, prim;
		int j, k;
		for (j = 0; j < size; j++) {
			temp = numbers[j];
			if (temp.signum() < 0) {
				temp = temp.negate();
				matrix[j][0] = 1;
			} else
				matrix[j][0] = 0;
			for (k = 1; k < PRIME_BASE; k++) {
				matrix[j][k] = 0;
				prim = bigPrimes[k];
				while (temp.mod(prim).compareTo(ZERO) == 0) {
					matrix[j][k]++; // = 1 - matrix[j][k];
					temp = temp.divide(prim);
				}
			}
			for (k = PRIME_BASE; k < size; k++)
				matrix[j][k] = 0;
		}
		return matrix;
	}

	private int[][] flattenMatrix(byte matrix[][]) {
		int m[][] = new int[matrix.length][matrix.length / flattenFactor];
		int j, k, n;
		int comparation;
		for (j = 0; j < matrix.length; j++)
			for (k = 0; k < matrix.length / flattenFactor; k++) {
				comparation = 1;
				m[j][k] = 0;
				for (n = 0; n < flattenFactor; n++) {
					if ((matrix[j][k * flattenFactor + n] & 1) > 0)
						m[j][k] += comparation;
					comparation *= 2;
				}
			}
		return m;
	}

	private int[][] buildIdentityMatrix(int size) {

		int matrix[][] = new int[size][size / flattenFactor];
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size / flattenFactor; j++)
				matrix[i][j] = 0;

		int j = -1;
		int comparison = 0;
		for (int i = 0; i < size; i++) {
			if (i % flattenFactor == 0) {
				j++;
				comparison = 1;
			} else
				comparison *= 2;
			matrix[i][j] = comparison;
		}
		return matrix;
	}

	private void divideByPower(int index, BigInteger p) {
		do
			fs[index] = fs[index].divide(p);
		while (fs[index].mod(p).equals(ZERO));
		if (fs[index].equals(ONE)) {
			logs[index] = 0;
			decomposedNumbers.add(s[index]);
			decomposedCount++;
		}
	}

	private long[] findQxPolynomials(long p, long n) {
		long k, x;

		long result[] = new long[2];
		result[0] = -1;
		result[1] = -1;

		if (p == 2) {
			// Special case. The prime is 2, and n will be either 0 or 1.
			result[0] = n;
			result[1] = -1;
			return result;
		}

		// Find a h value for which h^2-4n is not a quadratic residue mod p.
		long h = 3;
		do
			h += 2;
		while (legendreSymbol(h * h - 4 * n, p) != -1);

		// k is p/2 rounded up (p is always odd, of course).
		k = (p + 1) / 2;

		// Solves the equation of t congruent with n mod p.
		// Algorithm from D.M. Bressoud. Factorization and Primality Testing
		// 1989.
		x = v_(k, h, n, p);
		if (x < 0)
			x += p;
		x = (x * k) % p;
		result[0] = x;
		result[1] = (p - x);

		return result;
	}

	public static long v_(long j, long h, long n, long p) {
		long b[] = new long[64];
		long m = n;
		long v = h;
		long w = (h * h - 2 * m) % p;
		long x;
		int k, t;
		t = 0;
		while (j > 0) {
			b[++t] = j % 2;
			j /= 2;
		}
		for (k = t - 1; k >= 1; k--) {
			x = (v * w - h * m) % p;
			v = (v * v - 2 * m) % p;
			w = (w * w - 2 * n * m) % p;
			m = m * m % p;
			if (b[k] == 0)
				w = x;
			else {
				v = x;
				m = n * m % p;
			}
		}
		return v;
	}

	public static long legendreSymbol(long n, long p) {
		long count, temp;
		long legendre = 1;
		if (n == 0)
			return 0;
		if (n < 0) {
			n = -n;
			if (p % 4 == 3)
				legendre = -1;
		}
		do {
			count = 0;
			while (n % 2 == 0) {
				n = n / 2;
				count = 1 - count;
			}
			if ((count * (p * p - 1)) % 16 == 8)
				legendre = -legendre;
			if (((n - 1) * (p - 1)) % 8 == 4)
				legendre = -legendre;
			temp = n;
			n = p % n;
			p = temp;
		} while (n > 1);
		return legendre;
	}

	public static long modPowLong(long n, long p, long m) {
		if (p == 0)
			return 1;
		if (p % 2 == 1)
			return (n * modPowLong(n, p - 1, m)) % m;
		else {
			long result = modPowLong(n, p / 2, m);
			return (result * result) % m;
		}
	}

	public static BigInteger sqrt(BigInteger x) throws IllegalArgumentException {
		if (x.compareTo(BigInteger.ZERO) < 0) {
			throw new IllegalArgumentException("Negative argument.");
		}
		// square roots of 0 and 1 are trivial and
		// y == 0 will cause a divide-by-zero exception
		if (x.equals(BigInteger.ZERO) || x.equals(BigInteger.ONE)) {
			return x;
		} // end if
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger y;
		// starting with y = x / 2 avoids magnitude issues with x squared
		for (y = x.divide(two); y.compareTo(x.divide(y)) > 0; y = ((x.divide(y)).add(y)).divide(two))
			;
		if (x.compareTo(y.multiply(y)) == 0) {
			return y;
		} else {
			return y.add(BigInteger.ONE);
		}
	} // end bigIntSqRootCeil
}
