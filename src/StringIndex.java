public class StringIndex implements Comparable<StringIndex> {
	public String s;
	public int i;

	public StringIndex(String s, int i) {
		this.s = s;
		this.i = i;
	}

	@Override
	public int compareTo(StringIndex o) {
		return this.i - o.i;
	}
}