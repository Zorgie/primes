import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

public class Big {

	private static final int POLLARD_RHO = 0;
	private static final int FERMAT = 1;
	private static final int QUADRATIC_SIEVE = 2;

	private static final BigInteger ZERO = new BigInteger("0");
	private static final BigInteger ONE = new BigInteger("1");
	private static final BigInteger TWO = new BigInteger("2");
	public static final int TIMEOUT = 14700;
	public static long start;

	public static void main(String[] args) {
		ArrayList<StringIndex> input = new ArrayList<StringIndex>();
		start = System.currentTimeMillis();
		Scanner scanner = new Scanner(System.in);
		String line;
		int index = 0;
		while (index < 100 && !(line = scanner.nextLine()).equals("0")) {
			input.add(new StringIndex(line, index++));
		}
		HashMap<Integer, List<BigInteger>> mapping = new HashMap<Integer, List<BigInteger>>();
		Collections.sort(input, new SICmp());
		for (StringIndex s : input) {
			if (System.currentTimeMillis() > start + TIMEOUT)
				break;
			mapping.put(s.i, rhoPollardFactors(new BigInteger(s.s)));
		}
		scanner.close();

		for (int i = 0; i < input.size(); i++) {
			if (i > 0)
				System.out.println();
			List<BigInteger> facts = mapping.get(i);
			if (facts == null) {
				System.out.println("fail");
			} else {
				for (BigInteger bi : facts) {
					System.out.println(bi);
				}
			}
		}
	}

	public static BigInteger getFactor(BigInteger N) {
		return ONE;
	}

	public static List<BigInteger> fermatFactors(BigInteger N) {
		List<BigInteger> list = new LinkedList<BigInteger>();
		if (N.isProbablePrime(15)) {
			list.add(N);
			return list;
		}
		BigInteger oldN = new BigInteger("0");
		while (oldN != N) {
			oldN = N;
			if (!N.testBit(0)) {
				list.add(TWO);
				N = N.divide(TWO);
				continue;
			}
			BigInteger a = sqrt(N);
			BigInteger b2 = a.pow(2).subtract(N);
			while (!isSquare(b2) && a.compareTo(N) < 0) {
				if (System.currentTimeMillis() - TIMEOUT > start)
					return null;
				b2 = b2.add(a.multiply(TWO).add(ONE));
				a = a.add(ONE);
			}
			if (a.equals(N)) {
				return list;
			} else {
				BigInteger fact1 = a.add(sqrt(b2));// (a + sqrt(b2));
				BigInteger fact2 = a.subtract(sqrt(b2));// (a - sqrt(b2));
				if (fact1.equals(N) || fact2.equals(N)) {
					list.add(N);
					return list;
				}
				if (!fact1.equals(ONE)) {
					List<BigInteger> muchos1 = fermatFactors(fact1);
					if (muchos1 == null)
						return null;
					list.addAll(muchos1);
				}
				if (!fact2.equals(ONE)) {
					List<BigInteger> muchos2 = fermatFactors(fact2);
					if (muchos2 == null)
						return null;
					list.addAll(muchos2);
				}
				N = N.divide(fact1);
				N = N.divide(fact2);
			}
		}
		return list;
	}

	public static List<BigInteger> rhoPollardFactors(BigInteger N) {
		List<BigInteger> list = new LinkedList<BigInteger>();
		BigInteger xi = TWO;
		BigInteger x2i = TWO;
		int limit = 500000;
		int loops = 0;
		do {
			if (loops++ % 1000 == 0) {
				if (System.currentTimeMillis() - TIMEOUT > start)
					return null;
			}
			BigInteger xiPrime = xi.pow(2).add(ONE);
			BigInteger x2iPrime = (x2i.pow(2).add(ONE).pow(2).add(ONE));
			xi = xiPrime.mod(N);
			x2i = x2iPrime.mod(N);
			if (loops % 5 == 0) {
				BigInteger s = N.gcd(xi.subtract(x2i));
				if (!s.equals(ONE) && !s.equals(N)) {
					N = N.divide(s);
					list.add(s);
				}
			}
		} while (!N.isProbablePrime(50));
		list.add(N);
		return list;
	}

	private static boolean isSquare(BigInteger s) {
		BigInteger root = sqrt(s);
		BigInteger test = root.pow(2);
		return test.equals(s);
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