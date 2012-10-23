import java.util.Comparator;

public class SICmp implements Comparator<StringIndex> {
	@Override
	public int compare(StringIndex s1, StringIndex s2) {
		if (s1.s.length() > s2.s.length()) {
			return 1;
		} else if (s1.s.length() < s2.s.length()) {
			return -1;
		}
		return s1.s.compareTo(s2.s);
	}
}